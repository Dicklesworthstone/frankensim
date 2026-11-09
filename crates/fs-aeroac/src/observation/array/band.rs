//! Shared positive temporal factors for cross-location acoustic moments.

use super::*;
use fs_fft::harmonic_window::{
    HarmonicBandRef, HarmonicBandWorkspace, MAX_UNIFORM_HARMONIC_BAND_MODES,
    MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
type Components = [[Dd; 2]; 4];

// At most three velocity axes are admitted by PreparedHarmonicArray.
type Projected = [Modes; 4];

pub(super) fn accumulate<const D: usize>(
    left: &PreparedHarmonicField<D>,
    right: &PreparedHarmonicField<D>,
    window: HarmonicWindow,
    kind: FactorKind,
    moments: &mut CrossSums<D>,
    averaged_covariance: Option<&mut CrossSums<D>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, AeroacError> {
    if left.lines.len() > MAX_UNIFORM_HARMONIC_BAND_TONES
        || right.lines.len() > MAX_UNIFORM_HARMONIC_BAND_TONES
    {
        return Ok(false);
    }
    let Some(first) = left.lines.first() else {
        return Ok(false);
    };
    let angular = first.frequency.is_angular();
    let mut frequencies = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let mut count = 0;
    for line in left.lines.iter().chain(&right.lines) {
        poll(cancelled)?;
        if line.frequency.is_angular() != angular {
            return Ok(false);
        }
        let frequency = line.frequency.native();
        if !frequencies[..count].contains(&frequency) {
            if count == frequencies.len() {
                return Ok(false);
            }
            frequencies[count] = frequency;
            count += 1;
        }
    }
    if count < 3 {
        return Ok(false);
    }
    let frequencies = &mut frequencies[..count];
    if !core::ptr::eq(left, right) {
        // Canonical union order makes opposite receiver queries use the same
        // factor and reduction order, including different spectral supports.
        frequencies.sort_unstable_by(f64::total_cmp);
    }
    // A self pair retains the original authored order, exactly matching its
    // local observer rather than introducing a second summation permutation.
    let mut workspace = HarmonicBandWorkspace::new();
    let band = match (kind, angular) {
        (FactorKind::Deterministic, false) => {
            window.prepare_real_band_hz_in(frequencies, &mut workspace)
        }
        (FactorKind::Deterministic, true) => {
            window.prepare_real_band_rad_s_in(frequencies, &mut workspace)
        }
        (FactorKind::Proper, false) => window.prepare_band_hz_in(frequencies, &mut workspace),
        (FactorKind::Proper, true) => window.prepare_band_rad_s_in(frequencies, &mut workspace),
    }
    .map_err(window_error)?;
    poll(cancelled)?;
    let Some(band) = band else {
        return Ok(false);
    };
    let a = project(left, frequencies, &band, cancelled)?;
    let b = project(right, frequencies, &band, cancelled)?;
    let normalization = match kind {
        FactorKind::Deterministic => Dd::ONE,
        FactorKind::Proper => Dd::from_f64(0.5),
    };
    for mode in 0..band.mode_count() {
        poll(cancelled)?;
        add(
            moments,
            &band,
            mode,
            std::array::from_fn(|component| a[component][mode]),
            std::array::from_fn(|component| b[component][mode]),
            normalization,
        )?;
    }
    if let Some(covariance) = averaged_covariance {
        let (a, b) = match kind {
            FactorKind::Deterministic => (
                std::array::from_fn(|component| a[component][0]),
                std::array::from_fn(|component| b[component][0]),
            ),
            FactorKind::Proper => {
                // A proper analytic band removes one common carrier for
                // SECOND moments. Its constant mode is not the physical mean.
                // Average the original-clock fields first, preserving their
                // carrier and complex phase before the source outer product.
                let a = averaged_amplitudes(left, window, cancelled)?;
                let b = averaged_amplitudes(right, window, cancelled)?;
                (components(&a), components(&b))
            }
        };
        // Mode zero has weight one. The proper half enters before products,
        // so a finite covariance need not form a larger analytic norm first.
        add(covariance, &band, 0, a, b, normalization)?;
    }
    poll(cancelled)?;
    Ok(true)
}

fn project<const D: usize>(
    field: &PreparedHarmonicField<D>,
    frequencies: &[f64],
    band: &HarmonicBandRef<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Projected, AeroacError> {
    let mut groups = [0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (group, line) in groups.iter_mut().zip(&field.lines) {
        poll(cancelled)?;
        *group = frequencies
            .iter()
            .position(|&frequency| frequency == line.frequency.native())
            .ok_or_else(|| invalid("acoustic temporal union lost an original frequency"))?;
    }
    let mut projected = [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 4];
    for (component, output) in projected.iter_mut().enumerate().take(D + 1) {
        poll(cancelled)?;
        let mut fields = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
        for (&group, line) in groups.iter().zip(&field.lines) {
            poll(cancelled)?;
            fields[group] = if component == 0 {
                line.pressure.parts()
            } else {
                line.velocity[component - 1].parts()
            };
        }
        // A frequency absent from one receiver has zero field, not an
        // independent source or a fabricated phase.
        band.project_parts_dd(
            &fields[..frequencies.len()],
            &mut output[..band.mode_count()],
        )
        .map_err(window_error)?;
        poll(cancelled)?;
    }
    Ok(projected)
}

fn components<const D: usize>(line: &PreparedLine<D>) -> Components {
    let mut fields = [[Dd::ZERO; 2]; 4];
    fields[0] = line.pressure.parts();
    for (out, value) in fields[1..].iter_mut().zip(line.velocity) {
        *out = value.parts();
    }
    fields
}

fn add<const D: usize>(
    sums: &mut CrossSums<D>,
    band: &HarmonicBandRef<'_>,
    mode: usize,
    left: Components,
    right: Components,
    normalization: Dd,
) -> Result<(), AeroacError> {
    let symmetric = left == right;
    let product = |row: usize, column: usize| {
        // Match the local observer's upper-triangle multiplication order for
        // an exact self field, without symmetrizing a cross-location tensor.
        let (row, column) = if symmetric && row > column {
            (column, row)
        } else {
            (row, column)
        };
        band.mode_cross_parts_dd_scaled(mode, left[row], right[column], normalization)
            .map(|value| value[0])
            .map_err(window_error)
    };
    sums.pp = sums.pp + product(0, 0)?;
    for row in 0..D {
        sums.pu[row] = sums.pu[row] + product(0, row + 1)?;
        sums.up[row] = sums.up[row] + product(row + 1, 0)?;
        for column in 0..D {
            sums.uu[row][column] = sums.uu[row][column] + product(row + 1, column + 1)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
