//! Native carrier and relative-frequency phases for coherent causal taps.

use super::*;
use fs_fft::harmonic_window::{phase_rad_s_dd_from_clock_dd, phase_rad_s_from_clock_dd};

/// Cache the original carrier once per nonempty coherent bank and authored tap.
/// A delay advances the complex phasor under the exp(-i omega t) convention.
pub(super) fn carrier(omega: f64, delay: f64) -> Result<[Dd; 2], TransferError> {
    phase_rad_s_from_clock_dd(omega, Dd::from_f64(-delay)).map_err(window_error)
}

/// Apply the original gain after a separately retained relative phase. Nearby
/// tones share their carrier's phase rounding; independent bright trigonometric
/// values must not overwhelm the weak field that survives their interference.
pub(super) fn weight(
    reference: f64,
    omega: f64,
    tap: HarmonicPressureTap,
    [cr, ci]: [Dd; 2],
) -> Result<[Dd; 2], TransferError> {
    let offset = Dd::from_f64(omega) - Dd::from_f64(reference);
    let [dr, di] =
        phase_rad_s_dd_from_clock_dd(offset, Dd::from_f64(-tap.delay_s)).map_err(window_error)?;
    let weight = Dd::from_f64(tap.weight);
    let result = [(cr * dr - ci * di) * weight, (cr * di + ci * dr) * weight];
    if result.iter().any(|&value| !finite(value)) {
        Err(bad("pressure array delayed gain is unrepresentable"))
    } else {
        Ok(result)
    }
}
