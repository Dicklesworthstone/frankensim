//! Project coherent source fields before quadratic energy and its pullback.
//!
//! The shared temporal basis is unrotated here: immutable source groups
//! already carry the exact joint aperture/source-clock center phase. Both
//! complex parts are retained until the original spatial even/odd Gram acts.

use super::*;
use fs_fft::harmonic_window::{HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandWorkspace};

pub(in super::super) struct FieldAction {
    coefficients: Vec<Dd>,
    modes: usize,
    combined: Vec<FieldParts>,
    projected: Vec<FieldParts>,
    gradient: Vec<[Dd; 2]>,
}

impl FieldAction {
    /// Conservative extra logical payload and complete projected query work.
    /// Existing owner/Gram/source payload remains reserved by the caller.
    pub(in super::super) fn requirements(
        grid: AngularSpectrumGrid,
        frequencies: usize,
        source_groups: usize,
        controls: usize,
        modes: usize,
    ) -> Result<(usize, usize), AngularSpectrumError> {
        let n = mul(grid.width, grid.height)?;
        let work = mul(
            mul(n, modes)?,
            add(
                add(add(frequencies, source_groups)?, controls)?,
                add(add(grid.width, grid.height)?, 1)?,
            )?,
        )?;
        let bytes = add(
            add(
                mul(mul(frequencies, modes)?, core::mem::size_of::<Dd>())?,
                mul(
                    mul(n, add(frequencies, 1)?)?,
                    core::mem::size_of::<FieldParts>(),
                )?,
            )?,
            add(
                mul(controls, core::mem::size_of::<[Dd; 2]>())?,
                core::mem::size_of::<HarmonicBandWorkspace>()
                    + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES,
            )?,
        )?;
        Ok((work, bytes))
    }

    pub(in super::super) fn new(
        window: HarmonicWindow,
        frequencies: &[f64],
        spatial_modes: usize,
        controls: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Self>, AngularSpectrumError> {
        poll(cancelled)?;
        // The original groups already include the time-center rotation.
        let mut workspace = HarmonicBandWorkspace::new();
        let centered = HarmonicWindow::with_profile(0.0, window.duration_s(), window.profile())
            .map_err(|_| bad("optical control temporal interval is unrepresentable"))?;
        let Some(factor) = centered
            .prepare_real_band_rad_s_in(frequencies, &mut workspace)
            .map_err(|_| bad("optical control temporal factor is unrepresentable"))?
        else {
            return Ok(None);
        };
        let modes = factor.mode_count();
        let mut coefficients = reserve(mul(frequencies.len(), modes)?)?;
        for mode in 0..modes {
            poll(cancelled)?;
            for source in 0..frequencies.len() {
                coefficients.push(
                    factor
                        .centered_coefficient_dd(mode, source)
                        .map_err(|_| bad("optical control temporal coefficient is unavailable"))?,
                );
            }
        }
        let result = Self {
            coefficients,
            modes,
            combined: field_storage(mul(spatial_modes, frequencies.len())?, cancelled)?,
            projected: field_storage(spatial_modes, cancelled)?,
            gradient: filled(controls, [Dd::ZERO; 2])?,
        };
        poll(cancelled)?;
        Ok(Some(result))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate(
        &mut self,
        gram: &mut GramWorkspace,
        spatial: &Spatial<RealKernel>,
        groups: &[SourceGroup],
        weights: [Dd; 6],
        output: &mut [C64],
        include_energy: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        let n = self.projected.len();
        let f = self.coefficients.len() / self.modes;
        for (output, input) in self.combined.chunks_mut(128).zip(gram.residual.chunks(128)) {
            poll(cancelled)?;
            output.copy_from_slice(input);
        }
        self.gradient.fill([Dd::ZERO; 2]);
        let mut energy = Dd::ZERO;
        for temporal_mode in 0..self.modes {
            let coefficients = &self.coefficients[temporal_mode * f..(temporal_mode + 1) * f];
            for (mode, projected) in self.projected.iter_mut().enumerate() {
                poll(cancelled)?;
                *projected = [[Dd::ZERO; 2]; 6];
                for (frequency, &coefficient) in coefficients.iter().enumerate() {
                    let fields = &self.combined[frequency * n + mode];
                    for (out, field) in projected.iter_mut().zip(fields) {
                        out[0] = out[0] + coefficient * field[0];
                        out[1] = out[1] + coefficient * field[1];
                    }
                }
                if projected.iter().flatten().any(|&v| !finite_dd(v)) {
                    return Err(bad("optical control projected field is unrepresentable"));
                }
            }
            for (output, input) in gram.residual[..n]
                .chunks_mut(128)
                .zip(self.projected.chunks(128))
            {
                poll(cancelled)?;
                output.copy_from_slice(input);
            }
            // Odd temporal parity exchanges the spatial cosine/sine channel.
            // This is exactly Re(i*F), without discarding imaginary spatial
            // amplitudes or inferring transverse/optical-cycle semantics.
            let weight = Dd::from_f64((2 * temporal_mode + 1) as f64);
            let temporal = if temporal_mode % 2 == 0 {
                RealKernel([weight, Dd::ZERO])
            } else {
                RealKernel([Dd::ZERO, weight])
            };
            gram.real(spatial, &[temporal], 1, cancelled)?;
            for (field, covectors) in self.projected.iter().zip(&mut gram.covectors[..n]) {
                poll(cancelled)?;
                for ((value, covector), weight) in field.iter().zip(covectors).zip(weights) {
                    covector[0] = covector[0] * weight;
                    covector[1] = covector[1] * weight;
                    if covector.iter().any(|&v| !finite_dd(v)) {
                        return Err(bad("optical control field covector is unrepresentable"));
                    }
                    if include_energy {
                        energy = energy
                            + (value[0] * covector[0] + value[1] * covector[1]) * Dd::from_f64(0.5);
                    }
                }
            }
            // A source can span several frequencies. Project its complete
            // field BEFORE the pullback too; summing large per-frequency
            // gradients would reintroduce the cancellation just removed.
            let mut cursor = 0;
            for (control, derivative) in self.gradient.iter_mut().enumerate() {
                while cursor < groups.len() && groups[cursor].control.is_none() {
                    cursor += 1;
                }
                let start = cursor;
                while cursor < groups.len() && groups[cursor].control == Some(control) {
                    cursor += 1;
                }
                if start == cursor {
                    continue;
                }
                for mode in 0..n {
                    poll(cancelled)?;
                    let mut source = [[Dd::ZERO; 2]; 6];
                    for group in &groups[start..cursor] {
                        let coefficient = coefficients[group.frequency];
                        for (out, field) in source.iter_mut().zip(&group.fields[mode]) {
                            out[0] = out[0] + coefficient * field[0];
                            out[1] = out[1] + coefficient * field[1];
                        }
                    }
                    for ([fr, fi], [gr, gi]) in source.iter().zip(&gram.covectors[mode]) {
                        derivative[0] = derivative[0] + *fr * *gr + *fi * *gi;
                        derivative[1] = derivative[1] - *fi * *gr + *fr * *gi;
                    }
                }
            }
        }
        if !finite_dd(energy)
            || !energy.to_f64().is_finite()
            || energy.hi < 0.0
            || (energy.hi == 0.0 && energy.lo < 0.0)
        {
            return Err(bad("optical control energy density is unrepresentable"));
        }
        for (result, derivative) in output.iter_mut().zip(&self.gradient) {
            if derivative
                .iter()
                .any(|&v| !finite_dd(v) || !v.to_f64().is_finite())
            {
                return Err(bad("optical source gradient is unrepresentable"));
            }
            *result = C64::new(derivative[0].to_f64(), derivative[1].to_f64());
        }
        poll(cancelled)?;
        Ok(energy.to_f64())
    }
}

fn field_storage(
    count: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<FieldParts>, AngularSpectrumError> {
    let mut fields = reserve(count)?;
    for i in 0..count {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        fields.push([[Dd::ZERO; 2]; 6]);
    }
    Ok(fields)
}
