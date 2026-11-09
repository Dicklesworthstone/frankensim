//! Cached temporal field projection for repeated coherent arm observations.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, MAX_UNIFORM_HARMONIC_BAND_MODES,
    MAX_UNIFORM_HARMONIC_BAND_TONES, PreparedHarmonicBand,
};

type Source = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

// The owned factor's two Vec headers are part of Bank. Shared preparation
// and allocation-free read scratch are reserved once, separately from the
// retained factor arrays and mode-major projected spatial fields.
const FIXED_BYTES: usize = HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
    + 2 * core::mem::size_of::<Source>()
    + core::mem::size_of::<Modes>()
    + core::mem::size_of::<[f64; MAX_UNIFORM_HARMONIC_BAND_TONES]>();

pub(super) struct Plan {
    pub(super) modes: usize,
    pub(super) extra_work: usize,
    pub(super) extra_bytes: usize,
}

impl Plan {
    pub(super) fn new(
        grid: AngularSpectrumGrid,
        window: HarmonicWindow,
        frequencies: &[f64],
        proper: bool,
        reserved_modes: usize,
    ) -> Result<Option<Self>, AngularSpectrumError> {
        let modes = if proper {
            window.band_mode_count_rad_s(frequencies)
        } else {
            window.real_band_mode_count_rad_s(frequencies)
        }
        .map_err(|_| bad("aperture interferometer temporal admission failed"))?;
        let Some(modes) = modes else {
            return Ok(None);
        };
        let overflow = || bad("aperture interferometer projection bounds overflow");
        let n = grid.width.checked_mul(grid.height).ok_or_else(overflow)?;
        let axes = grid.width.checked_add(grid.height).ok_or_else(overflow)?;
        let f = frequencies.len();
        let old_work = f
            .checked_mul(n)
            .and_then(|v| axes.checked_add(f).and_then(|passes| v.checked_mul(passes)))
            .ok_or_else(overflow)?;
        let new_work = modes
            .checked_mul(n)
            .and_then(|v| {
                axes.checked_add(f)
                    .and_then(|p| p.checked_add(1))
                    .and_then(|passes| v.checked_mul(passes))
            })
            .and_then(|v| {
                if proper {
                    n.checked_mul(f)
                        .and_then(|rotations| v.checked_add(rotations))
                } else {
                    Some(v)
                }
            })
            .ok_or_else(overflow)?;
        let cache_bytes = f
            .checked_mul(modes + 2)
            .and_then(|v| v.checked_mul(core::mem::size_of::<Dd>()))
            .ok_or_else(overflow)?;
        let more_modes = modes.saturating_sub(reserved_modes);
        let more_fields = n
            .checked_mul(more_modes)
            .and_then(|v| v.checked_mul(core::mem::size_of::<FieldParts>()))
            .ok_or_else(overflow)?;
        let extra_bytes = cache_bytes
            .checked_add(more_fields)
            .and_then(|v| v.checked_add(if reserved_modes == 0 { FIXED_BYTES } else { 0 }))
            .ok_or_else(overflow)?;
        Ok(Some(Self {
            modes,
            extra_work: new_work.saturating_sub(old_work),
            extra_bytes,
        }))
    }
}

pub(super) struct Positive {
    factor: PreparedHarmonicBand,
    proper: bool,
}

impl Positive {
    pub(super) fn new(
        window: HarmonicWindow,
        frequencies: &[f64],
        proper: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        let mut negative = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
        let factor = if proper {
            if frequencies.len() > negative.len() {
                return Err(bad("aperture interferometer temporal source count changed"));
            }
            for (out, &frequency) in negative.iter_mut().zip(frequencies) {
                *out = -frequency;
            }
            window.prepare_band_rad_s(&negative[..frequencies.len()])
        } else {
            // Real source groups already carry the joint spatial/clock phase.
            // Their unrotated coefficients preserve both complex spatial parts.
            HarmonicWindow::with_profile(0.0, window.duration_s(), window.profile())
                .map_err(|_| bad("aperture interferometer centered interval is unrepresentable"))?
                .prepare_real_band_rad_s(frequencies)
        }
        .map_err(|_| bad("aperture interferometer temporal factor is unrepresentable"))?
        .ok_or_else(|| bad("aperture interferometer temporal factor changed after admission"))?;
        poll(cancelled)?;
        Ok(Self { factor, proper })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn accumulate(
        &self,
        gram: &mut GramWorkspace,
        real: Option<&Spatial<RealKernel>>,
        proper: Option<&Spatial<ProperKernel>>,
        points: usize,
        projected: &mut [FieldParts],
        moments: &mut [Dd; 21],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        let modes = self.factor.mode_count();
        let f = self.factor.source_count();
        if projected.len() < points * modes || gram.residual.len() < points * f {
            return Err(bad(
                "aperture interferometer projected storage is inconsistent",
            ));
        }
        for point in 0..points {
            poll(cancelled)?;
            if self.proper {
                for component in 0..6 {
                    let mut source: Source = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
                    for (frequency, value) in source[..f].iter_mut().enumerate() {
                        *value = gram.residual[frequency * points + point][component];
                    }
                    let mut output: Modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
                    self.factor
                        .project_parts_dd(&source[..f], &mut output[..modes])
                        .map_err(|_| {
                            bad("aperture interferometer temporal field is unrepresentable")
                        })?;
                    for mode in 0..modes {
                        projected[mode * points + point][component] = output[mode];
                    }
                }
            } else {
                for mode in 0..modes {
                    let output = &mut projected[mode * points + point];
                    *output = [[Dd::ZERO; 2]; 6];
                    for source in 0..f {
                        let coefficient = self
                            .factor
                            .centered_coefficient_dd(mode, source)
                            .map_err(|_| {
                                bad("aperture interferometer coefficient is unavailable")
                            })?;
                        for (out, input) in output
                            .iter_mut()
                            .zip(gram.residual[source * points + point])
                        {
                            out[0] = out[0] + coefficient * input[0];
                            out[1] = out[1] + coefficient * input[1];
                        }
                    }
                }
            }
        }
        for tile in projected[..points * modes].chunks(128) {
            poll(cancelled)?;
            if tile.iter().flatten().flatten().any(|&v| !finite_dd(v)) {
                return Err(bad(
                    "aperture interferometer temporal field is unrepresentable",
                ));
            }
        }
        for mode in 0..modes {
            let fields = &projected[mode * points..(mode + 1) * points];
            for (output, input) in gram.residual[..points]
                .chunks_mut(128)
                .zip(fields.chunks(128))
            {
                poll(cancelled)?;
                output.copy_from_slice(input);
            }
            if self.proper {
                gram.proper(
                    proper.ok_or_else(|| bad("missing proper interferometer cache"))?,
                    &[ProperKernel::temporal_mode(mode)],
                    1,
                    cancelled,
                )?;
            } else {
                gram.real(
                    real.ok_or_else(|| bad("missing real interferometer cache"))?,
                    &[RealKernel::temporal_mode(mode)],
                    1,
                    cancelled,
                )?;
            }
            for (field, covector) in fields.iter().zip(&gram.covectors[..points]) {
                poll(cancelled)?;
                let mut pair = 0;
                for i in 0..6 {
                    for j in i..6 {
                        moments[pair] = moments[pair]
                            + field[i][0] * covector[j][0]
                            + field[i][1] * covector[j][1];
                        pair += 1;
                    }
                }
            }
        }
        poll(cancelled)
    }
}

#[cfg(test)]
mod tests;
