//! Coherent linear channel readouts retaining their fast absolute carrier.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_READOUT_QUERY_SCRATCH_BYTES, HarmonicReadoutRef, HarmonicReadoutWorkspace,
    MAX_HARMONIC_READOUT_TERMS, MAX_HARMONIC_READOUT_TONES,
};

pub(super) type Fields = [[Dd; 2]; MAX_HARMONIC_READOUT_TONES];

pub(super) fn eligible(count: usize) -> bool {
    (2..=MAX_HARMONIC_READOUT_TONES).contains(&count)
}

/// The original linear contractions or a complete bounded readout preparation
/// and all requested complex means, whichever is larger. This work is added
/// to second-moment work; their fixed scratch is used sequentially.
pub(super) fn admission(count: usize, outputs: usize) -> Result<(usize, usize), TransferError> {
    let original = count
        .checked_mul(outputs)
        .and_then(|p| p.checked_mul(2))
        .ok_or_else(overflow)?;
    if outputs == 0 || !eligible(count) {
        return Ok((original, 0));
    }
    let preparation = count
        .checked_mul(MAX_HARMONIC_READOUT_TERMS + 4)
        .and_then(|p| p.checked_add(16 * MAX_HARMONIC_READOUT_TERMS + 4))
        .ok_or_else(overflow)?;
    let projection = count
        .checked_add(1)
        .and_then(|p| p.checked_mul(4 + 2 * MAX_HARMONIC_READOUT_TERMS))
        .and_then(|p| p.checked_mul(outputs))
        .and_then(|p| p.checked_add(preparation))
        .ok_or_else(overflow)?;
    let scratch = core::mem::size_of::<HarmonicReadoutWorkspace>()
        + core::mem::size_of::<[f64; MAX_HARMONIC_READOUT_TONES]>()
        + core::mem::size_of::<Fields>()
        + HARMONIC_READOUT_QUERY_SCRATCH_BYTES;
    Ok((original.max(projection), scratch))
}

pub(super) fn prepare<'a>(
    bank: &PressureBank,
    window: HarmonicWindow,
    workspace: &'a mut HarmonicReadoutWorkspace,
) -> Result<Option<HarmonicReadoutRef<'a>>, TransferError> {
    if !eligible(bank.lines.len()) {
        return Ok(None);
    }
    let mut frequencies = [0.0; MAX_HARMONIC_READOUT_TONES];
    for (out, line) in frequencies.iter_mut().zip(&bank.lines) {
        *out = line.omega;
    }
    window
        .prepare_high_carrier_readout_rad_s_in(&frequencies[..bank.lines.len()], workspace)
        .map_err(window_error)
}

/// False leaves all means unchanged for the original absolute-clock kernels.
pub(super) fn accumulate(
    bank: &PressureBank,
    window: HarmonicWindow,
    integrated: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    if integrated.is_empty() || !eligible(bank.lines.len()) {
        return Ok(false);
    }
    poll(cancelled)?;
    let mut workspace = HarmonicReadoutWorkspace::new();
    let Some(factor) = prepare(bank, window, &mut workspace)? else {
        return Ok(false);
    };
    let mut fields: Fields = [[Dd::ZERO; 2]; MAX_HARMONIC_READOUT_TONES];
    for (channel, mean) in integrated.iter_mut().enumerate() {
        poll(cancelled)?;
        for (field, line) in fields.iter_mut().zip(&bank.lines) {
            *field = line.amplitudes[channel];
        }
        *mean = factor
            .mean_parts_dd(&fields[..bank.lines.len()])
            .map_err(window_error)?;
    }
    poll(cancelled)?;
    Ok(true)
}
