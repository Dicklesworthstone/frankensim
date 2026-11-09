//! Positive temporal modes on the original continuous spatial Fourier fields.
//!
//! The real path receives the original joint clock/aperture-centered phasors.
//! Both complex parts survive the frequency projection: real-field temporal
//! parity composes with spatial parity only when the final product is taken.
//! Proper factors retain their original complex spatial difference kernels.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandWorkspace,
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Source = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

// Shared preparation, one source/returned projection and its internal source
// scratch, signed frequencies, and the retained-vector header. Scalar locals,
// allocator overhead and compiler call-stack layout are not counted.
const FIXED_BYTES: usize = core::mem::size_of::<HarmonicBandWorkspace>()
    + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
    + 2 * core::mem::size_of::<Source>()
    + core::mem::size_of::<Modes>()
    + core::mem::size_of::<[f64; MAX_UNIFORM_HARMONIC_BAND_TONES]>()
    + core::mem::size_of::<Vec<FieldParts>>();

pub(super) fn admit(
    window: HarmonicWindow,
    frequencies: &[f64],
    points: usize,
    proper: bool,
    work_headroom: &mut usize,
    workspace_headroom: usize,
) -> Result<Option<usize>, AngularSpectrumError> {
    let modes = if proper {
        window.band_mode_count_rad_s(frequencies)
    } else {
        window.real_band_mode_count_rad_s(frequencies)
    }
    .map_err(|_| bad("aperture temporal projection admission failed"))?;
    let Some(modes) = modes else {
        return Ok(None);
    };
    let Some(bytes) = points
        .checked_mul(modes)
        .and_then(|v| v.checked_mul(core::mem::size_of::<FieldParts>()))
        .and_then(|v| v.checked_add(FIXED_BYTES))
    else {
        return Ok(None);
    };
    let Some(work) = points
        .checked_add(frequencies.len())
        .and_then(|v| v.checked_mul(points))
        .and_then(|v| v.checked_mul(modes))
    else {
        return Ok(None);
    };
    let Some(pair_work) = points
        .checked_mul(frequencies.len())
        .and_then(|v| v.checked_mul(v))
    else {
        return Ok(None);
    };
    let extra = work.saturating_sub(pair_work);
    if bytes > workspace_headroom || extra > *work_headroom {
        return Ok(None);
    }
    *work_headroom -= extra;
    Ok(Some(modes))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn integrate<K: Kernel, A: Accumulator>(
    spatial: &Spatial<K>,
    window: HarmonicWindow,
    frequencies: &[f64],
    fields: &[Vec<FieldParts>],
    modes: usize,
    sum: &mut A,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let points = spatial.grid.width * spatial.grid.height;
    let projected = project::<K>(window, frequencies, fields, points, modes, cancelled)?;
    for mode in 0..modes {
        let temporal = K::temporal_mode(mode);
        for i in 0..points {
            poll(cancelled)?;
            for j in 0..points {
                if j % 32 == 0 {
                    poll(cancelled)?;
                }
                sum.add(
                    spatial.pair(i, j).compose(temporal),
                    &projected[mode * points + i],
                    &projected[mode * points + j],
                );
            }
        }
    }
    poll(cancelled)
}

/// Mode-major original spatial fields. First spatial moments share the same
/// linear projection before signed force/torque contractions.
pub(super) fn project<K: Kernel>(
    window: HarmonicWindow,
    frequencies: &[f64],
    fields: &[Vec<FieldParts>],
    points: usize,
    modes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<FieldParts>, AngularSpectrumError> {
    poll(cancelled)?;
    let mut workspace = HarmonicBandWorkspace::new();
    let mut negative = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let factor = if K::CENTER_FIELDS {
        // Center rotations have already been applied jointly with space. Only
        // translation-independent coefficients are used; avoiding a second
        // source-clock rotation preserves translated spatiotemporal nulls.
        HarmonicWindow::with_profile(0.0, window.duration_s(), window.profile())
            .map_err(|_| bad("aperture centered temporal window is unrepresentable"))?
            .prepare_real_band_rad_s_in(frequencies, &mut workspace)
    } else {
        if frequencies.len() > negative.len() {
            return Err(bad(
                "aperture temporal source count changed after admission",
            ));
        }
        for (out, &frequency) in negative.iter_mut().zip(frequencies) {
            *out = -frequency;
        }
        // EM uses exp(+i omega t). The common carrier cancels only here,
        // within this one authored proper-complex source bank.
        window.prepare_band_rad_s_in(&negative[..frequencies.len()], &mut workspace)
    }
    .map_err(|_| bad("aperture temporal factor is unrepresentable"))?
    .ok_or_else(|| bad("aperture temporal factor changed after admission"))?;
    if modes != factor.mode_count() || fields.len() != factor.source_count() {
        return Err(bad("aperture temporal factor shape is inconsistent"));
    }
    let count = points
        .checked_mul(modes)
        .ok_or_else(|| bad("aperture temporal field count overflow"))?;
    let mut projected = Vec::new();
    projected
        .try_reserve_exact(count)
        .map_err(|_| bad("aperture temporal field allocation failed"))?;
    for i in 0..count {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        projected.push([[Dd::ZERO; 2]; 6]);
    }
    if fields.iter().any(|field| field.len() != points) {
        return Err(bad("aperture temporal field shape is inconsistent"));
    }
    for point in 0..points {
        poll(cancelled)?;
        if K::CENTER_FIELDS {
            for mode in 0..modes {
                let out = &mut projected[mode * points + point];
                for (source, field) in fields.iter().enumerate() {
                    let coefficient = factor
                        .centered_coefficient_dd(mode, source)
                        .map_err(|_| bad("aperture temporal coefficient is unrepresentable"))?;
                    for (output, input) in out.iter_mut().zip(field[point]) {
                        for part in 0..2 {
                            output[part] = output[part] + coefficient * input[part];
                        }
                    }
                }
            }
        } else {
            for component in 0..6 {
                let mut source: Source = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
                for (value, field) in source.iter_mut().zip(fields) {
                    *value = field[point][component];
                }
                let mut output: Modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
                factor
                    .project_parts_dd(&source[..frequencies.len()], &mut output[..modes])
                    .map_err(|_| bad("aperture temporal field is unrepresentable"))?;
                for mode in 0..modes {
                    projected[mode * points + point][component] = output[mode];
                }
            }
        }
    }
    for tile in projected.chunks(128) {
        poll(cancelled)?;
        if tile
            .iter()
            .flatten()
            .flatten()
            .any(|&value| !finite_dd(value))
        {
            return Err(bad("aperture temporal field is unrepresentable"));
        }
    }
    poll(cancelled)?;
    Ok(projected)
}

#[cfg(test)]
mod tests;
