//! Exterior stress and field momentum on one original temporal field factor.
//!
//! All inputs use the shared exp(-i omega t) convention. Spatial cosine/sine
//! and carrier Gram representations are linear in those same amplitudes, so
//! temporal projection commutes with the original continuous cell integral.
//! Uniform and triangular observations use their actual normalized sensitivity;
//! a weighted momentum rate is not an unweighted endpoint change divided by time.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandRef, HarmonicBandWorkspace,
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

struct Projected {
    value: [Modes; 6],
    derivative: [Modes; 6],
}

// The carrier path retains four projections. Include return temporaries,
// source/derivative arrays, temporal preparation and the projection owner's
// own fixed input scratch conservatively. No heap is added to an observation.
pub(in super::super) const WORKSPACE_BYTES: usize = core::mem::size_of::<HarmonicBandWorkspace>()
    + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
    + 6 * core::mem::size_of::<Projected>()
    + 4 * core::mem::size_of::<[[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES]>()
    + core::mem::size_of::<[Dd; MAX_UNIFORM_HARMONIC_BAND_TONES]>();

struct Temporal<'a> {
    factor: HarmonicBandRef<'a>,
    derivative_frequencies: [Dd; MAX_UNIFORM_HARMONIC_BAND_TONES],
    proper: bool,
    stress: bool,
    momentum: bool,
    epsilon: Dd,
    mu: Dd,
}

impl Temporal<'_> {
    fn project(
        &self,
        field: impl Fn(usize) -> [[Dd; 2]; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Projected, AngularSpectrumError> {
        let mut projected = Projected {
            value: [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 6],
            derivative: [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 6],
        };
        let count = self.factor.source_count();
        let modes = self.factor.mode_count();
        for component in 0..6 {
            poll(cancelled)?;
            let mut source = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
            let mut derivative = source;
            for i in 0..count {
                source[i] = field(i)[component];
                if self.momentum {
                    let [re, im] = source[i];
                    let omega = self.derivative_frequencies[i];
                    derivative[i] = [im * omega, -re * omega];
                }
            }
            self.factor
                .project_parts_dd(&source[..count], &mut projected.value[component][..modes])
                .map_err(|_| bad("exterior temporal field projection is unrepresentable"))?;
            if self.momentum {
                self.factor
                    .project_parts_dd(
                        &derivative[..count],
                        &mut projected.derivative[component][..modes],
                    )
                    .map_err(|_| bad("exterior temporal derivative is unrepresentable"))?;
            }
        }
        poll(cancelled)?;
        Ok(projected)
    }

    fn cross(&self, left: &Modes, right: &Modes, scale: Dd) -> Result<Dd, AngularSpectrumError> {
        let scale = scale
            * if self.proper {
                Dd::from_f64(0.5)
            } else {
                Dd::ONE
            };
        let mut value = Dd::ZERO;
        for mode in 0..self.factor.mode_count() {
            value = value
                + self
                    .factor
                    .mode_cross_parts_dd_scaled(mode, left[mode], right[mode], scale)
                    .map_err(|_| bad("exterior temporal cross moment is unrepresentable"))?[0];
        }
        if !value.hi.is_finite() || !value.lo.is_finite() {
            return Err(bad("exterior temporal cross moment overflows"));
        }
        Ok(value)
    }

    fn moments(
        &self,
        left: &Projected,
        right: &Projected,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Dd; 6], [Dd; 6]), AngularSpectrumError> {
        let mut stress = [Dd::ZERO; 6];
        let mut momentum = [Dd::ZERO; 6];
        if self.stress {
            for (c, (u, v)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
                .into_iter()
                .enumerate()
            {
                poll(cancelled)?;
                stress[c] = self.cross(&left.value[u], &right.value[v], self.epsilon)?
                    + self.cross(&left.value[u + 3], &right.value[v + 3], self.mu)?;
            }
        }
        if self.momentum {
            for c in 0..3 {
                poll(cancelled)?;
                let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                for (e, h, sign) in [(u, v + 3, Dd::ONE), (v, u + 3, -Dd::ONE)] {
                    momentum[c] = momentum[c]
                        + sign * self.cross(&left.value[e], &right.value[h], Dd::ONE)?;
                    momentum[c + 3] = momentum[c + 3]
                        + sign * self.cross(&left.derivative[e], &right.value[h], Dd::ONE)?
                        + sign * self.cross(&left.value[e], &right.derivative[h], Dd::ONE)?;
                }
            }
        }
        poll(cancelled)?;
        Ok((stress, momentum))
    }

    fn cell(
        &self,
        grid: AngularSpectrumGrid,
        fields: &[Vec<[C64; 6]>],
        carrier: Option<&cell::CarrierFields>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Dd; 6], [Dd; 6]), AngularSpectrumError> {
        let mut stress = [Dd::ZERO; 6];
        let mut momentum = [Dd::ZERO; 6];
        let area = Dd::from_f64(grid.width as f64 * grid.spacing_m[0])
            * Dd::from_f64(grid.height as f64 * grid.spacing_m[1]);
        for mode in 0..grid.width * grid.height {
            poll(cancelled)?;
            if let Some(carrier) = carrier {
                // The original separable spatial Gram pass supplies source,
                // cosine-weighted and sine-weighted complex coefficients.
                // Its second product uses i*A, which must be projected in
                // its own real temporal quadrature before contraction.
                let source =
                    self.project(|i| carrier.at(i, mode).map(|c| c.parts()[0]), cancelled)?;
                let cosine =
                    self.project(|i| carrier.at(i, mode).map(|c| c.parts()[1]), cancelled)?;
                let source_i = self.project(
                    |i| carrier.at(i, mode).map(|c| rotate_i(c.parts()[0])),
                    cancelled,
                )?;
                let sine_i = self.project(
                    |i| carrier.at(i, mode).map(|c| rotate_i(c.parts()[2])),
                    cancelled,
                )?;
                for (a, b) in [(&source, &cosine), (&source_i, &sine_i)] {
                    let (dyad, flow) = self.moments(a, b, cancelled)?;
                    add(&mut stress, dyad, area);
                    add(&mut momentum, flow, area);
                }
                continue;
            }
            let opposite = opposite_mode(grid, mode);
            if opposite.is_some_and(|other| other < mode) {
                continue;
            }
            let constant = opposite == Some(mode);
            let partner = opposite.filter(|_| !constant);
            let weight = area * Dd::from_f64(if constant { 1.0 } else { 0.5 });
            for quadrature in 0..if constant { 1 } else { 2 } {
                let field = self.project(
                    |i| {
                        spatial_quadratures(fields[i][mode], partner.map(|other| fields[i][other]))
                            [quadrature]
                    },
                    cancelled,
                )?;
                let (dyad, flow) = self.moments(&field, &field, cancelled)?;
                add(&mut stress, dyad, weight);
                add(&mut momentum, flow, weight);
            }
        }
        if stress
            .iter()
            .chain(&momentum)
            .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
        {
            return Err(bad("exterior temporal cell integral overflows"));
        }
        poll(cancelled)?;
        Ok((stress, momentum))
    }
}

fn rotate_i([re, im]: [Dd; 2]) -> [Dd; 2] {
    [-im, re]
}

fn add(target: &mut [Dd; 6], value: [Dd; 6], weight: Dd) {
    for c in 0..6 {
        target[c] = target[c] + value[c] * weight;
    }
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn try_accumulate(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    proper: bool,
    epsilon: f64,
    mu: f64,
    mut stress: Option<&mut [[Dd; 6]]>,
    mut momentum: Option<&mut [[Dd; 6]]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<([Dd; 6], [Dd; 6])>, AngularSpectrumError> {
    poll(cancelled)?;
    let mut scratch = HarmonicBandWorkspace::new();
    let factor = if proper {
        window.prepare_band_rad_s_in(&bank.frequencies, &mut scratch)
    } else {
        window.prepare_real_band_rad_s_in(&bank.frequencies, &mut scratch)
    }
    .map_err(|_| bad("exterior temporal factor is unrepresentable"))?;
    let Some(factor) = factor else {
        return Ok(None);
    };
    let mut derivative_frequencies = [Dd::ZERO; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let reference = if proper {
        Dd::from_f64(bank.frequencies[0])
    } else {
        Dd::ZERO
    };
    for (out, &omega) in derivative_frequencies.iter_mut().zip(&bank.frequencies) {
        *out = Dd::from_f64(omega) - reference;
    }
    let temporal = Temporal {
        factor,
        derivative_frequencies,
        proper,
        stress: stress.is_some(),
        momentum: momentum.is_some(),
        epsilon: Dd::from_f64(epsilon),
        mu: Dd::from_f64(mu),
    };
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(bank.frequency_count())
        .map_err(|_| bad("exterior temporal field allocation failed"))?;
    for group in 0..bank.frequency_count() {
        poll(cancelled)?;
        fields.push(bank.frequency_modal_field(group, side, distance, cancelled)?);
    }
    let carrier = bank.carrier_fields(&fields, proper, cancelled)?;
    let integrated = temporal.cell(bank.anchor.grid(), &fields, carrier.as_ref(), cancelled)?;
    drop(carrier);
    for (group, field) in fields.iter_mut().enumerate() {
        bank.synthesize_frequency_field(group, field, cancelled)?;
    }
    for point in 0..bank.anchor.modes.len() {
        poll(cancelled)?;
        let projected = temporal.project(
            |i| fields[i][point].map(|z| [Dd::from_f64(z.re), -Dd::from_f64(z.im)]),
            cancelled,
        )?;
        let (dyad, flow) = temporal.moments(&projected, &projected, cancelled)?;
        if let Some(sums) = stress.as_deref_mut() {
            add(&mut sums[point], dyad, Dd::ONE);
        }
        if let Some(sums) = momentum.as_deref_mut() {
            add(&mut sums[point], flow, Dd::ONE);
        }
    }
    poll(cancelled)?;
    Ok(Some(integrated))
}

#[cfg(test)]
mod tests;
