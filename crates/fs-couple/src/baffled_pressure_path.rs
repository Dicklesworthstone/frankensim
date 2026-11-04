//! Causal path convolution shared by lossless and material-loaded Rayleigh kernels.
//! This leaf owns interpolation/composition only, never a propagation law.
use fs_math::det;
use std::f64::consts::{PI, TAU};

pub(super) const HALF: usize = 32;
pub(super) const TAPS: usize = 2 * HALF + 1;

pub(super) fn window() -> [f64; TAPS] {
    let mut w = [0.0; TAPS];
    for (i, value) in w.iter_mut().enumerate().take(TAPS - 1).skip(1) {
        let phase = TAU * i as f64 / (TAPS - 1) as f64;
        *value = 0.42 - 0.5 * det::cos(phase) + 0.08 * det::cos(2.0 * phase);
    }
    w
}

pub(super) fn taps(fraction: f64, window: &[f64; TAPS]) -> [f64; TAPS] {
    let mut h = [0.0; TAPS];
    if fraction == 0.0 {
        h[HALF] = 1.0;
        return h;
    }
    let sine = det::sin(PI * fraction.min(1.0 - fraction));
    for (i, value) in h.iter_mut().enumerate().take(TAPS - 1).skip(1) {
        let offset = i as i32 - HALF as i32;
        let numerator = if offset % 2 == 0 { -sine } else { sine };
        *value = window[i] * numerator / (PI * (f64::from(offset) - fraction));
    }
    let sum: f64 = h.iter().sum();
    for value in &mut h {
        *value /= sum;
    }
    h
}

/// Convolve one source-derived causal response with its OWN residual travel
/// delay. The response is neither normalized nor replaced by its magnitude.
/// `delay_samples` excludes the response's own delay and includes any common
/// explicit padding selected by the caller; every emitted sample is causal.
pub(super) fn compose(
    response: &[f64],
    delay_samples: f64,
    maximum_frames: usize,
    remaining_terms: &mut usize,
) -> Result<(usize, Vec<f64>), String> {
    if response.is_empty()
        || response.iter().any(|x| !x.is_finite())
        || !delay_samples.is_finite()
        || delay_samples < 0.0
        || delay_samples >= maximum_frames as f64
    {
        return Err("a Rayleigh path needs a finite causal response and residual delay".into());
    }
    let count = response
        .len()
        .checked_add(TAPS - 1)
        .ok_or("Rayleigh path response extent overflow")?;
    let first = delay_samples.floor() as usize;
    if first
        .checked_add(count)
        .is_none_or(|end| end > maximum_frames)
    {
        return Err("material propagation exceeds the Rayleigh history budget".into());
    }
    let work = response
        .len()
        .checked_mul(TAPS)
        .ok_or("Rayleigh path convolution work overflow")?;
    if work > *remaining_terms {
        return Err("material propagation exceeds its declared convolution work".into());
    }
    let fractional = taps(delay_samples - first as f64, &window());
    let mut combined = vec![0.0; count];
    for (i, &source) in response.iter().enumerate() {
        for (j, &weight) in fractional.iter().enumerate() {
            combined[i + j] += source * weight;
        }
    }
    if combined.iter().any(|x| !x.is_finite()) {
        return Err("material propagation path convolution overflow".into());
    }
    *remaining_terms -= work;
    Ok((first, combined))
}

#[cfg(test)]
#[path = "baffled_pressure_path_tests.rs"]
mod tests;
