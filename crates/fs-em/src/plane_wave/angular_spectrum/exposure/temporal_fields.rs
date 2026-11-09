//! Original exterior E/H fields projected before normalized-window products.
//!
//! Uniform and triangular sensitivity use their own positive temporal factor.
//! Real fields retain the full carrier and all sum tones; authored proper
//! factors use half of the analytic cross moment. Spatial cosine/sine and
//! carrier Gram fields remain the original continuous-cell representation.
//! These bounded optional projections require separately reserved headroom;
//! unsupported windows and legacy budgets keep the original pair law.

use super::*;
pub(super) use fs_fft::harmonic_window::HarmonicBandWorkspace as Workspace;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandRef, MAX_UNIFORM_HARMONIC_BAND_MODES,
    MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
type Projected = [Modes; 6];

// Six vector projections cover retained values and return temporaries; two
// source arrays cover the caller and shared projection scratch. Logical payload
// excludes scalar locals, compiler temporaries and allocator/call-stack overhead.
pub(super) const WORKSPACE_BYTES: usize = core::mem::size_of::<Workspace>()
    + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
    + 6 * core::mem::size_of::<Projected>()
    + 2 * core::mem::size_of::<[[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES]>();

fn eligible(frequencies: usize, proper: bool) -> bool {
    (3..=MAX_UNIFORM_HARMONIC_BAND_TONES / if proper { 1 } else { 2 }).contains(&frequencies)
}

/// Optional reservation, decided before field reconstruction or query allocation.
/// Four six-component vector passes cover point moments and the carrier cell's
/// four four-component projections, including their modal cross contractions.
/// Source-count/span refusal never changes another source bank's admission.
pub(super) fn admit(
    points: usize,
    frequencies: usize,
    proper: bool,
    work_headroom: usize,
    workspace_headroom: usize,
) -> Option<usize> {
    if !eligible(frequencies, proper) || workspace_headroom < WORKSPACE_BYTES {
        return None;
    }
    let work = points
        .checked_mul(frequencies.checked_add(1)?)?
        .checked_mul(4 * MAX_UNIFORM_HARMONIC_BAND_MODES)?;
    (work != 0 && work <= work_headroom).then_some(work)
}

pub(super) fn admit_banks(
    points: usize,
    bank_count: usize,
    banks: impl IntoIterator<Item = (usize, bool)>,
    mut work_headroom: usize,
    workspace_headroom: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<bool>, usize, usize), AngularSpectrumError> {
    let mut flags = Vec::new();
    let mut terms = 0;
    let Some(flag_bytes) = bank_count.checked_mul(core::mem::size_of::<bool>()) else {
        return Ok((flags, 0, 0));
    };
    let Some(headroom) = workspace_headroom.checked_sub(flag_bytes) else {
        return Ok((flags, 0, 0));
    };
    for (index, (frequencies, proper)) in banks.into_iter().enumerate() {
        poll(cancelled)?;
        if index >= bank_count {
            return Err(bad("exterior projection bank count is inconsistent"));
        }
        let work = admit(points, frequencies, proper, work_headroom, headroom);
        if work.is_some() && terms == 0 {
            flags
                .try_reserve_exact(bank_count)
                .map_err(|_| bad("exterior projection admission allocation failed"))?;
            flags.resize(index, false);
        }
        if terms != 0 || work.is_some() {
            flags.push(work.is_some());
        }
        if let Some(work) = work {
            work_headroom -= work;
            terms += work;
        }
    }
    poll(cancelled)?;
    let bytes = if terms == 0 {
        0
    } else {
        WORKSPACE_BYTES + flag_bytes
    };
    Ok((flags, terms, bytes))
}

pub(super) struct PositiveFields<'w> {
    factor: HarmonicBandRef<'w>,
    proper: bool,
}

impl<'w> PositiveFields<'w> {
    pub(super) fn prepare(
        window: HarmonicWindow,
        frequencies: &[f64],
        proper: bool,
        workspace: &'w mut Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Self>, AngularSpectrumError> {
        poll(cancelled)?;
        if !eligible(frequencies.len(), proper) {
            return Ok(None);
        }
        let factor = if proper {
            window.prepare_band_rad_s_in(frequencies, workspace)
        } else {
            window.prepare_real_band_rad_s_in(frequencies, workspace)
        }
        .map_err(|_| bad("exterior field temporal factor is unrepresentable"))?;
        poll(cancelled)?;
        Ok(factor.map(|factor| Self { factor, proper }))
    }

    fn project(
        &self,
        components: &[usize],
        field: impl Fn(usize) -> [[Dd; 2]; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Projected, AngularSpectrumError> {
        let mut projected = [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 6];
        let count = self.factor.source_count();
        for &component in components {
            poll(cancelled)?;
            let mut source = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
            for (i, value) in source[..count].iter_mut().enumerate() {
                *value = field(i)[component];
            }
            self.factor
                .project_parts_dd(
                    &source[..count],
                    &mut projected[component][..self.factor.mode_count()],
                )
                .map_err(|_| bad("exterior field temporal projection is unrepresentable"))?;
        }
        poll(cancelled)?;
        Ok(projected)
    }

    fn cross(&self, left: &Modes, right: &Modes) -> Result<Dd, AngularSpectrumError> {
        let scale = if self.proper {
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
                    .map_err(|_| bad("exterior field temporal product is unrepresentable"))?[0];
        }
        Ok(value)
    }

    fn normal_flux(&self, left: &Projected, right: &Projected) -> Result<Dd, AngularSpectrumError> {
        Ok(self.cross(&left[0], &right[4])? - self.cross(&left[1], &right[3])?)
    }

    pub(super) fn point(
        &self,
        field: impl Fn(usize) -> [C64; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 9], AngularSpectrumError> {
        let projected = self.project(
            &[0, 1, 2, 3, 4, 5],
            |i| field(i).map(|z| [Dd::from_f64(z.re), -Dd::from_f64(z.im)]),
            cancelled,
        )?;
        let mut moments = [Dd::ZERO; 9];
        for c in 0..6 {
            poll(cancelled)?;
            moments[c] = self.cross(&projected[c], &projected[c])?;
        }
        for c in 0..3 {
            let (u, v) = ((c + 1) % 3, (c + 2) % 3);
            moments[c + 6] = self.cross(&projected[u], &projected[v + 3])?
                - self.cross(&projected[v], &projected[u + 3])?;
        }
        if moments.iter().any(|v| !v.to_f64().is_finite()) {
            return Err(bad("exterior field moment is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(moments)
    }

    pub(super) fn cell_normal_power(
        &self,
        grid: AngularSpectrumGrid,
        fields: &[Vec<[C64; 6]>],
        carrier: Option<&cell::CarrierFields>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AngularSpectrumError> {
        let area = Dd::from_f64(grid.width as f64 * grid.spacing_m[0])
            * Dd::from_f64(grid.height as f64 * grid.spacing_m[1]);
        let mut power = Dd::ZERO;
        for mode in 0..grid.width * grid.height {
            poll(cancelled)?;
            if let Some(carrier) = carrier {
                // Stream the two original spatial Gram products. Multiplication
                // by i is performed before the real temporal projection.
                for (part, rotate) in [(1, false), (2, true)] {
                    let source = self.project(
                        &[0, 1, 3, 4],
                        |i| {
                            carrier
                                .at(i, mode)
                                .map(|c| rotate_parts(c.parts()[0], rotate))
                        },
                        cancelled,
                    )?;
                    let weighted = self.project(
                        &[0, 1, 3, 4],
                        |i| {
                            carrier
                                .at(i, mode)
                                .map(|c| rotate_parts(c.parts()[part], rotate))
                        },
                        cancelled,
                    )?;
                    power = power + self.normal_flux(&source, &weighted)?;
                }
                continue;
            }
            let opposite = opposite_mode(grid, mode);
            if opposite.is_some_and(|other| other < mode) {
                continue;
            }
            let constant = opposite == Some(mode);
            let partner = opposite.filter(|_| !constant);
            let weight = Dd::from_f64(if constant { 1.0 } else { 0.5 });
            let mut flux = Dd::ZERO;
            for quadrature in 0..if constant { 1 } else { 2 } {
                let field = self.project(
                    &[0, 1, 3, 4],
                    |i| {
                        spatial_quadratures(fields[i][mode], partner.map(|other| fields[i][other]))
                            [quadrature]
                    },
                    cancelled,
                )?;
                flux = flux + self.normal_flux(&field, &field)?;
            }
            power = power + flux * weight;
        }
        // Apply physical area once after the compensated spatial sum, as the
        // original pair path does. Per-mode area scaling can erase a finite
        // total when individual weighted modes are subnormal.
        let power = power * area;
        if !power.hi.is_finite() || !power.lo.is_finite() {
            return Err(bad("exterior field cell power is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(power)
    }
}

fn rotate_parts([real, imaginary]: [Dd; 2], rotate: bool) -> [Dd; 2] {
    if rotate {
        [-imaginary, real]
    } else {
        [real, imaginary]
    }
}

#[cfg(test)]
mod tests;
