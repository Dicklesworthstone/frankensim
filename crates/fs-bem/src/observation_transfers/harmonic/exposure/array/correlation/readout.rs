//! Absolute-carrier means for selected signed-lag readout covariance.

use super::super::readout as linear;
use super::*;
use fs_fft::harmonic_window::{HarmonicReadoutWorkspace, MAX_HARMONIC_READOUT_TONES};

pub(super) fn admission(
    count: usize,
    pairs: usize,
    first_channels: usize,
) -> Result<(usize, usize), TransferError> {
    let original = count
        .checked_mul(pairs)
        .and_then(|p| p.checked_mul(4))
        .ok_or_else(overflow)?;
    let outputs = first_channels.checked_add(pairs).ok_or_else(overflow)?;
    let (products, scratch) = linear::admission(count, outputs)?;
    let cache = if scratch == 0 {
        0
    } else {
        core::mem::size_of::<[[Dd; 2]; MAX_PRESSURE_ARRAY_CHANNELS]>()
            + core::mem::size_of::<[bool; MAX_PRESSURE_ARRAY_CHANNELS]>()
    };
    Ok((original.max(products), scratch + cache))
}

pub(super) fn accumulate(
    bank: &PressureBank,
    lagged: &CorrelationBank,
    pairs: &[HarmonicPressureCorrelationPair],
    window: HarmonicWindow,
    integrated: &mut [[[Dd; 2]; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    if pairs.is_empty() || !linear::eligible(bank.lines.len()) {
        return Ok(false);
    }
    poll(cancelled)?;
    let mut workspace = HarmonicReadoutWorkspace::new();
    let Some(factor) = linear::prepare(bank, window, &mut workspace)? else {
        return Ok(false);
    };
    let mut fields: linear::Fields = [[Dd::ZERO; 2]; MAX_HARMONIC_READOUT_TONES];
    let mut first = [[Dd::ZERO; 2]; MAX_PRESSURE_ARRAY_CHANNELS];
    let mut cached = [false; MAX_PRESSURE_ARRAY_CHANNELS];
    for (p, (pair, mean)) in pairs.iter().zip(integrated).enumerate() {
        poll(cancelled)?;
        if !cached[pair.first] {
            for (field, line) in fields.iter_mut().zip(&bank.lines) {
                *field = line.amplitudes[pair.first];
            }
            first[pair.first] = factor
                .mean_parts_dd(&fields[..bank.lines.len()])
                .map_err(window_error)?;
            cached[pair.first] = true;
        }
        for (line, field) in fields[..bank.lines.len()].iter_mut().enumerate() {
            *field = lagged.shifted[line * pairs.len() + p];
        }
        *mean = [
            first[pair.first],
            factor
                .mean_parts_dd(&fields[..bank.lines.len()])
                .map_err(window_error)?,
        ];
    }
    poll(cancelled)?;
    Ok(true)
}
