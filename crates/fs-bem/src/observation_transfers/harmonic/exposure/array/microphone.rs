//! Exact-frequency, calibrated pressure/velocity microphone data.
//!
//! A row is a settled linear response on the original exp(-i omega t) clock.
//! No interpolation, capsule dynamics, causality, noise or loading is inferred.

use super::*;

/// Complex sensitivity to the original Cartesian pressure/velocity field.
///
/// The scalar output is calibrated Pa-equivalent:
/// m = pressure_gain*p + sum_c particle_velocity_gain_pa_s_m[c]*u[c].
/// Pressure sensitivity is dimensionless; velocity sensitivity is Pa s/m.
/// Components use the ORIGINAL receiver geometry's Cartesian frame. The
/// caller must express every row in the same fixed authored orientation.
/// Complex sensitivities retain relative phase and a frequency-dependent local
/// first-order polar response; replacing actual u by p/(rho*c) changes the
/// measured quantity.
///
/// These are response data, not a claim of a causal capsule model. A finite
/// table establishes neither a causal filter nor behavior between its lines.
/// Electrical voltage/sensitivity must be explicitly calibrated to these units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicMicrophoneResponse {
    pub pressure_gain: C64,
    pub particle_velocity_gain_pa_s_m: [C64; 3],
}

/// One complete native-frequency calibration row, in physical receiver order.
///
/// Rows must be strictly increasing and match ALL distinct retained positive
/// angular frequencies exactly, including silent fields, unused receivers and
/// dark channels. There is no nearest-frequency matching or extrapolation.
/// Receiver identities/positions are those of the original BEM owner.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicMicrophoneCalibrationLine {
    pub angular_frequency_rad_s: f64,
    pub responses: Vec<HarmonicMicrophoneResponse>,
}

impl<'a> PreparedHarmonicPressureArray<'a> {
    /// Apply declared complex microphone response before coherent delayed taps.
    ///
    /// Mean and real/proper factors retain their original joint laws. The
    /// prepared array owns a copy of the complete calibration and exposes it
    /// beside moments/correlations. Original source and calibration uncertainty
    /// remain unresolved; pressure-equivalent output is not a new field point.
    pub fn new_with_calibration(
        exposure: &PreparedHarmonicExposure<'a>,
        calibration: &[HarmonicMicrophoneCalibrationLine],
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
    ) -> Result<Self, TransferError> {
        Self::new_with_calibration_cancellable(exposure, calibration, channels, limits, || false)
    }

    /// Atomic preparation with the existing complete-array resource bounds.
    pub fn new_with_calibration_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        calibration: &[HarmonicMicrophoneCalibrationLine],
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        Self::prepare(
            PressureInput::from_exposure(exposure),
            channels,
            None,
            Some(calibration),
            limits,
            cancelled,
        )
    }
}

// Shape admission precedes coefficient scans and allocation. The existing
// exposure cap bounds every row even before the complete work count is known.
pub(super) fn admit_shape(
    calibration: &[HarmonicMicrophoneCalibrationLine],
    receivers: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    admit(calibration.len(), MAX_EXPOSURE_LINES)?;
    for row in calibration {
        poll(cancelled)?;
        if row.responses.len() != receivers {
            return Err(bad("calibration needs every original physical receiver"));
        }
    }
    Ok(())
}

pub(super) fn extra_work(
    rows: usize,
    receivers: usize,
    lines: usize,
) -> Result<usize, TransferError> {
    // Sixteen real DD products per four-component complex projection. Bound
    // validation/conversion of all eight response scalars on every row too.
    receivers
        .checked_mul(lines.checked_mul(16).ok_or_else(overflow)?)
        .and_then(|v| v.checked_add(rows.checked_mul(receivers)?.checked_mul(16)?))
        .and_then(|v| v.checked_add(rows))
        .ok_or_else(overflow)
}

pub(super) fn extra_bytes(rows: usize, receivers: usize) -> Result<usize, TransferError> {
    // Retained data and per-frequency complex projection scratch. The common
    // projected scalar receiver array is counted by the pressure-array owner.
    rows.checked_mul(core::mem::size_of::<HarmonicMicrophoneCalibrationLine>())
        .and_then(|v| {
            v.checked_add(
                rows.checked_mul(receivers)?
                    .checked_mul(core::mem::size_of::<HarmonicMicrophoneResponse>())?,
            )
        })
        .and_then(|v| v.checked_add(receivers.checked_mul(core::mem::size_of::<[Dd; 8]>())?))
        .ok_or_else(overflow)
}

pub(super) fn copy_calibration(
    input: &[HarmonicMicrophoneCalibrationLine],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<HarmonicMicrophoneCalibrationLine>, TransferError> {
    let mut output = storage(input.len())?;
    for (index, row) in input.iter().enumerate() {
        poll(cancelled)?;
        if !row.angular_frequency_rad_s.is_finite()
            || row.angular_frequency_rad_s <= 0.0
            || (index != 0
                && row.angular_frequency_rad_s <= input[index - 1].angular_frequency_rad_s)
        {
            return Err(bad(
                "calibration frequencies must be finite, positive and strictly increasing",
            ));
        }
        let mut responses = storage(row.responses.len())?;
        for response in &row.responses {
            poll(cancelled)?;
            if core::iter::once(&response.pressure_gain)
                .chain(&response.particle_velocity_gain_pa_s_m)
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err(bad(
                    "microphone pressure/velocity sensitivities must be finite",
                ));
            }
            responses.push(*response);
        }
        output.push(HarmonicMicrophoneCalibrationLine {
            angular_frequency_rad_s: row.angular_frequency_rad_s,
            responses,
        });
    }
    Ok(output)
}

pub(super) fn admit_frequencies(
    calibration: &[HarmonicMicrophoneCalibrationLine],
    frequencies: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    if calibration.len() != frequencies.len() {
        return Err(bad(
            "calibration must cover exactly every retained physical frequency",
        ));
    }
    for (row, &omega) in calibration.iter().zip(frequencies) {
        poll(cancelled)?;
        if row.angular_frequency_rad_s.to_bits() != omega.to_bits() {
            return Err(bad(
                "microphone calibration differs from native angular frequency",
            ));
        }
    }
    Ok(())
}

// Convert once per unique frequency, shared across all mean/factor banks.
// No response lookup or coefficient conversion occurs during later sampling.
pub(super) fn coefficients(
    row: &HarmonicMicrophoneCalibrationLine,
    output: &mut [[Dd; 8]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    for (response, weights) in row.responses.iter().zip(output) {
        poll(cancelled)?;
        for (c, value) in core::iter::once(&response.pressure_gain)
            .chain(&response.particle_velocity_gain_pa_s_m)
            .enumerate()
        {
            weights[2 * c] = Dd::from_f64(value.re);
            weights[2 * c + 1] = Dd::from_f64(value.im);
        }
    }
    Ok(())
}

pub(super) fn project(
    input: &Line,
    projections: &[[Dd; 8]],
    selected: &[bool; crate::near_field::MAX_RECEIVERS],
    output: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    for (receiver, ((field, weights), signal)) in
        input.fields.iter().zip(projections).zip(output).enumerate()
    {
        poll(cancelled)?;
        if !selected[receiver] {
            continue;
        }
        *signal = [Dd::ZERO; 2];
        for c in 0..4 {
            let [ar, ai] = component(field, c);
            let (wr, wi) = (weights[2 * c], weights[2 * c + 1]);
            signal[0] = signal[0] + ar * wr - ai * wi;
            signal[1] = signal[1] + ar * wi + ai * wr;
        }
        if signal.iter().any(|&value| !finite(value)) {
            return Err(bad("calibrated microphone projection is unrepresentable"));
        }
    }
    Ok(())
}
