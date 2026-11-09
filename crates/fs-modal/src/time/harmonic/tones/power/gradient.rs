//! Coherent work/heat covectors through the original paired QR factors.

use super::{MAX_COORDINATES, ModalAcousticState, ModalAcousticTimeError, Sum, TonePower, invalid};
use crate::time::harmonic::tones::{MAX_HARMONIC_TONES, PreparedModalTones};

/// Physical derivatives of the two independent finite-interval power integrals.
/// Units are joules per unit of the input named by the containing field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalTonePowerDerivative {
    pub input_work: f64,
    pub viscous_dissipation: f64,
}

/// One original constant-envelope source, in the interval-start carrier basis.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalTonePowerSourceGradient {
    pub cosine_force: ModalTonePowerDerivative,
    pub sine_force: ModalTonePowerDerivative,
    /// Common phase in C*cos(Omega*t+phi) + S*sin(Omega*t+phi), at phi=0.
    pub initial_phase: ModalTonePowerDerivative,
}

/// One coherent mode's power gradients, retaining every original source.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalTonePowerGradient {
    pub initial_displacement: ModalTonePowerDerivative,
    pub initial_velocity: ModalTonePowerDerivative,
    sources: [ModalTonePowerSourceGradient; MAX_HARMONIC_TONES],
    count: usize,
}

impl ModalTonePowerGradient {
    /// Original caller order. Coincident frequencies remain separate controls.
    pub fn source_gradients(&self) -> &[ModalTonePowerSourceGradient] {
        &self.sources[..self.count]
    }
}

/// Both affine-envelope knots use the SAME interval-start carrier basis.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalAffineTonePowerSourceGradient {
    pub cosine_start: ModalTonePowerDerivative,
    pub sine_start: ModalTonePowerDerivative,
    pub cosine_end: ModalTonePowerDerivative,
    pub sine_end: ModalTonePowerDerivative,
    /// One common additional phase applied to both endpoint pairs.
    pub initial_phase: ModalTonePowerDerivative,
}

/// Coherent power gradients for one original mode with affine excitation.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalAffineTonePowerGradient {
    pub initial_displacement: ModalTonePowerDerivative,
    pub initial_velocity: ModalTonePowerDerivative,
    sources: [ModalAffineTonePowerSourceGradient; MAX_HARMONIC_TONES],
    count: usize,
}

impl ModalAffineTonePowerGradient {
    /// Original caller order, including individually controlled duplicates.
    pub fn source_gradients(&self) -> &[ModalAffineTonePowerSourceGradient] {
        &self.sources[..self.count]
    }
}

/// Explicit additional admission for a single-mode power-gradient query.
/// The existing bank owns every numerical coefficient; no new factor is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModalTonePowerGradientBudget {
    /// Conservative fixed-array payload bound, including result and nested
    /// source/transpose scratch. Excludes caller storage and stack-frame overhead.
    pub maximum_scratch_bytes: usize,
    /// Fixed algebraic loop score, 64*n^2 + 256*input_tones, with n=2+2K
    /// or n=2+4K for K distinct frequencies. This is not a timing guarantee.
    pub maximum_query_work: usize,
}

impl ModalTonePowerGradientBudget {
    #[must_use]
    pub const fn realtime_reference() -> Self {
        Self {
            maximum_scratch_bytes: 64 * 1024,
            maximum_query_work: 2_000_000,
        }
    }
}

/// Borrowed, admitted work/heat transpose of the original coherent tone bank.
///
/// State, all source amplitudes/knots and independent source phase are supported.
/// Mechanics k/c, source frequencies, duration and basis/mass normalization stay
/// fixed. Warm queries use O(n^2) factor contractions, allocate nothing and call
/// no transcendental functions. No single-tone heat summation or subtraction of
/// endpoint energies is used. Finite-arithmetic gradients are not enclosures.
#[derive(Debug)]
pub struct PreparedModalTonePowerGradient<'a> {
    bank: &'a PreparedModalTones,
}

impl PreparedModalTones {
    /// Admit a borrowed power-gradient query surface for the original factors.
    /// No numerical preparation, coefficient allocation or state change occurs.
    /// Cancellation/refusal publishes no partial result and preserves the bank.
    pub fn prepare_power_gradients(
        &self,
        budget: ModalTonePowerGradientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedModalTonePowerGradient<'_>, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        admit_layout(self.unique.len(), self.tone_count(), self.affine, budget)?;
        poll(&mut cancelled)?;
        Ok(PreparedModalTonePowerGradient { bank: self })
    }
}

impl PreparedModalTonePowerGradient<'_> {
    /// Coherent signed-work/positive-heat gradients for one constant-tone mode.
    /// Source rows have exactly the original tone count. Each independent phase
    /// derivative rotates its own prescribed quadratures before coalescing.
    pub fn gradients_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalTonePowerGradient, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if self.bank.affine {
            return Err(invalid(
                "affine power gradients require both envelope endpoints",
            ));
        }
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let (c, s) = self.bank.collect_forces(cosine, sine)?;
        let k = self.bank.unique.len();
        let gradient = mode
            .power
            .input_gradients(state, &c[..k], &s[..k], None, &mut cancelled)?;
        let mut sources = [ModalTonePowerSourceGradient::default(); MAX_HARMONIC_TONES];
        for (source, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let i = 2 + 2 * unique;
            sources[source] = ModalTonePowerSourceGradient {
                cosine_force: gradient[i],
                sine_force: gradient[i + 1],
                initial_phase: phase(gradient[i], gradient[i + 1], cosine[source], sine[source])?,
            };
        }
        poll(&mut cancelled)?;
        Ok(ModalTonePowerGradient {
            initial_displacement: gradient[0],
            initial_velocity: gradient[1],
            sources,
            count: self.bank.tone_count(),
        })
    }

    /// Both affine endpoint pairs are in the interval-start carrier basis.
    /// Every work/heat cross term remains, including initial transients and
    /// interference between different or coincident frequencies. A reverse
    /// trajectory adds these covectors to the motion pullback, sums adjacent
    /// knot contributions and transposes accepted carrier-basis rotations.
    #[allow(clippy::too_many_arguments)]
    pub fn gradients_mode_affine(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalAffineTonePowerGradient, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if !self.bank.affine {
            return Err(invalid(
                "affine power gradients require affine tone preparation",
            ));
        }
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let [c, s, dc, ds] =
            self.bank
                .collect_affine(cosine_start, sine_start, cosine_end, sine_end)?;
        let k = self.bank.unique.len();
        let gradient = mode.power.input_gradients(
            state,
            &c[..k],
            &s[..k],
            Some((&dc[..k], &ds[..k])),
            &mut cancelled,
        )?;
        let mut sources = [ModalAffineTonePowerSourceGradient::default(); MAX_HARMONIC_TONES];
        for (source, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let i = 2 + 4 * unique;
            let start_phase = phase(
                gradient[i],
                gradient[i + 1],
                cosine_start[source],
                sine_start[source],
            )?;
            let end_phase = phase(
                gradient[i + 2],
                gradient[i + 3],
                cosine_end[source],
                sine_end[source],
            )?;
            sources[source] = ModalAffineTonePowerSourceGradient {
                cosine_start: gradient[i],
                sine_start: gradient[i + 1],
                cosine_end: gradient[i + 2],
                sine_end: gradient[i + 3],
                initial_phase: finite(ModalTonePowerDerivative {
                    input_work: start_phase.input_work + end_phase.input_work,
                    viscous_dissipation: start_phase.viscous_dissipation
                        + end_phase.viscous_dissipation,
                })?,
            };
        }
        poll(&mut cancelled)?;
        Ok(ModalAffineTonePowerGradient {
            initial_displacement: gradient[0],
            initial_velocity: gradient[1],
            sources,
            count: self.bank.tone_count(),
        })
    }
}

impl TonePower {
    pub(super) fn input_gradients(
        &self,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        deltas: Option<(&[f64], &[f64])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[ModalTonePowerDerivative; MAX_COORDINATES], ModalAcousticTimeError> {
        let x = self.coordinates(state, cosine, sine, deltas)?;
        let n = 2 + (if self.affine { 4 } else { 2 }) * cosine.len();
        let mut work = [Sum::default(); MAX_COORDINATES];
        let mut heat = work;
        // Work = scale * (R*x)^T(T*x), heat = c*h * ||R*x||^2.
        // Apply the original factors before forming products, preserving small
        // coherent residuals. Transpose each physical input column within the
        // row before accumulating, rather than subtracting full large gradients.
        for (velocity, force) in self
            .heat_root
            .chunks_exact(n)
            .zip(self.force_rows.chunks_exact(n))
        {
            poll(cancelled)?;
            let (mut v, mut f) = (Sum::default(), Sum::default());
            for ((&vr, &fr), &coordinate) in velocity.iter().zip(force).zip(&x[..n]) {
                v.add(vr * coordinate);
                f.add(fr * coordinate);
            }
            let (wv, wf, hv) = (
                self.work_scale * v.value(),
                self.work_scale * f.value(),
                self.heat_scale * v.value(),
            );
            for j in 0..n {
                let dv = self.physical_column(velocity, j);
                let df = self.physical_column(force, j);
                work[j].add((self.work_scale * dv) * wf);
                work[j].add(wv * (self.work_scale * df));
                let term = (self.heat_scale * dv) * hv;
                heat[j].add(term);
                heat[j].add(term);
            }
        }
        let mut result = [ModalTonePowerDerivative::default(); MAX_COORDINATES];
        for i in 0..n {
            result[i] = finite(ModalTonePowerDerivative {
                input_work: work[i].value(),
                viscous_dissipation: heat[i].value(),
            })?;
        }
        poll(cancelled)?;
        Ok(result)
    }

    // Columns of the ORIGINAL source chart [v,(sum C-kq)/kappa,C/kappa,
    // S/kappa,DeltaC/kappa,DeltaS/kappa,...], transposed to [q,v,C0,S0,C1,S1].
    fn physical_column(&self, row: &[f64], column: usize) -> f64 {
        match column {
            0 => -((row[1] * self.omega) * (self.omega * self.h)) / self.scale,
            1 => row[0],
            _ => {
                let width = if self.affine { 4 } else { 2 };
                let component = (column - 2) % width;
                let mut sum = Sum::default();
                if component == 0 {
                    sum.add(row[1]);
                }
                sum.add(row[column]);
                if self.affine && component < 2 {
                    sum.add(-row[column + 2]);
                }
                (self.h / self.scale) * sum.value()
            }
        }
    }
}

fn phase(
    c: ModalTonePowerDerivative,
    s: ModalTonePowerDerivative,
    cosine: f64,
    sine: f64,
) -> Result<ModalTonePowerDerivative, ModalAcousticTimeError> {
    finite(ModalTonePowerDerivative {
        input_work: sine.mul_add(c.input_work, -cosine * s.input_work),
        viscous_dissipation: sine.mul_add(c.viscous_dissipation, -cosine * s.viscous_dissipation),
    })
}

pub(super) fn finite(
    value: ModalTonePowerDerivative,
) -> Result<ModalTonePowerDerivative, ModalAcousticTimeError> {
    if value.input_work.is_finite() && value.viscous_dissipation.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "coherent tone power gradients must be representable",
        ))
    }
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ModalAcousticTimeError> {
    if cancelled() {
        Err(invalid("tone power gradient query cancelled"))
    } else {
        Ok(())
    }
}

// Shared by direct input covectors and source-clock companion queries. The
// latter use the affine layout even though their nominal bank is constant.
pub(super) fn admit_layout(
    distinct: usize,
    input_tones: usize,
    affine: bool,
    budget: ModalTonePowerGradientBudget,
) -> Result<(), ModalAcousticTimeError> {
    let bytes = core::mem::size_of::<ModalTonePowerGradient>()
        + core::mem::size_of::<ModalAffineTonePowerGradient>()
        + core::mem::size_of::<PreparedModalTonePowerGradient<'_>>()
        + 4 * core::mem::size_of::<[ModalTonePowerDerivative; MAX_COORDINATES]>()
        + 4 * core::mem::size_of::<[f64; MAX_COORDINATES]>()
        + 4 * core::mem::size_of::<[Sum; MAX_COORDINATES]>()
        + 4 * core::mem::size_of::<[[f64; MAX_HARMONIC_TONES]; 4]>();
    cap(
        "tone power gradient scratch bytes",
        bytes,
        budget.maximum_scratch_bytes,
    )?;
    let n = 2 + (if affine { 4 } else { 2 }) * distinct;
    cap(
        "tone power gradient query work score",
        64 * n * n + 256 * input_tones,
        budget.maximum_query_work,
    )
}

fn cap(what: &'static str, value: usize, limit: usize) -> Result<(), ModalAcousticTimeError> {
    if value > limit {
        Err(ModalAcousticTimeError::BudgetExceeded {
            what,
            value: value as f64,
            limit: limit as f64,
        })
    } else {
        Ok(())
    }
}
