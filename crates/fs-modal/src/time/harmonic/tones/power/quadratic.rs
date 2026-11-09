//! Complementary quadratic companions for original affine source clocks.

use super::frequency::ModalTonePowerFrequencyGradient;
use super::gradient::{ModalTonePowerDerivative, ModalTonePowerGradientBudget, finite, poll};
use super::{
    ModalAcousticState, ModalAcousticTimeError, Sum, TonePower, invalid, scaling_with_degree,
};
use crate::time::harmonic::tones::{MAX_HARMONIC_TONES, ModalToneBudget, PreparedModalTones};

const MAX_QUADRATIC_COORDINATES: usize = 2 + 6 * MAX_HARMONIC_TONES;

pub(super) fn admit(
    bank: &PreparedModalTones,
    preparation: ModalToneBudget,
    query: ModalTonePowerGradientBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ModalAcousticTimeError> {
    let k = bank.unique.len();
    let n = (2 + 6 * k) as u128;
    let square = n * n;
    // Two retained paired factors/mode; one extra dense multiplication/flow
    // workspace, scalar convolution/endpoint moments, and small row buffers.
    // The original borrowed bank, allocator metadata and caller inputs exclude.
    let storage = bank.mode_count() as u128
        * (core::mem::size_of::<TonePower>() as u128 + 16 * square)
        + 32 * square
        + 8 * (6 * n + 32 * k as u128)
        + 512;
    cap(
        "quadratic power companion storage bytes",
        storage,
        preparation.maximum_storage_bytes,
    )?;
    cap(
        "quadratic power companion step work score",
        bank.mode_count() as u128 * (16 * square + 128 * bank.tone_count() as u128),
        preparation.maximum_step_work,
    )?;
    let mut work = 0u128;
    for mode in &bank.modes {
        poll(cancelled)?;
        let motion = &mode.motion[0];
        let (_, _, _, _, doubles) = scaling_with_degree(
            motion.angular_frequency_rad_s,
            motion.transition.drag,
            &bank.unique,
            bank.duration_s,
            2,
        )?;
        work += 192 * square
            + doubles as u128 * (8 * n * square + 80 * square)
            + k as u128 * (12_000 + 512 * doubles as u128);
        cap(
            "quadratic power companion preparation work score",
            work,
            preparation.maximum_preparation_work,
        )?;
    }
    // All arrays are fixed at the hard source ceiling; include source
    // coalescing and returned derivatives in the nested warm payload bound.
    let scratch = 4 * core::mem::size_of::<[f64; MAX_QUADRATIC_COORDINATES]>()
        + 4 * core::mem::size_of::<[Sum; MAX_HARMONIC_TONES]>()
        + 4 * core::mem::size_of::<[[f64; MAX_HARMONIC_TONES]; 4]>()
        + 2 * core::mem::size_of::<ModalTonePowerFrequencyGradient>();
    cap(
        "quadratic power frequency scratch bytes",
        scratch as u128,
        query.maximum_scratch_bytes,
    )?;
    cap(
        "quadratic power frequency query work score",
        32 * square + 128 * n * bank.tone_count() as u128,
        query.maximum_query_work,
    )?;
    poll(cancelled)
}

fn cap(what: &'static str, value: u128, limit: usize) -> Result<(), ModalAcousticTimeError> {
    if value > limit as u128 {
        Err(ModalAcousticTimeError::BudgetExceeded {
            what,
            value: value as f64,
            limit: limit as f64,
        })
    } else {
        Ok(())
    }
}

impl TonePower {
    pub(super) fn try_new_quadratic(
        omega: f64,
        drag: f64,
        frequencies: &[f64],
        h: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        Self::prepare(omega, drag, frequencies, h, 2)
    }

    // The two native clock columns are u*(1-u) and u*u. The exact source
    // direction is h*[u*(1-u)*(S0,-C0) + u*u*(S1,-C1)]. Contract the ORIGINAL
    // source endpoint pairs per factor row; no individual endpoint delta and
    // no subtraction of nearly equal full covectors is formed.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn affine_clock_derivatives(
        &self,
        state: ModalAcousticState,
        nominal: [&[f64]; 4],
        endpoints: [&[f64]; 4],
        source_to_unique: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[ModalTonePowerDerivative; MAX_HARMONIC_TONES], ModalAcousticTimeError> {
        if self.envelope_degree != 2 {
            return Err(invalid(
                "affine source clocks require quadratic power companions",
            ));
        }
        let [cosine, sine, dc, ds] = nominal;
        let [c0, s0, c1, s1] = endpoints;
        let n = 2 + 6 * cosine.len();
        let mut x = [0.0; MAX_QUADRATIC_COORDINATES];
        self.coordinates_into(state, cosine, sine, Some((dc, ds)), &mut x)?;
        let mut work = [Sum::default(); MAX_HARMONIC_TONES];
        let mut heat = work;
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
            for (source, &unique) in source_to_unique.iter().enumerate() {
                let j = 4 + 6 * unique;
                let clock = |row: &[f64]| {
                    let mut sum = Sum::default();
                    // Form the physical column scales before large endpoint
                    // amplitudes, preserving admitted source cancellation.
                    sum.add(((self.h / self.scale) * row[j]) * s0[source]);
                    sum.add(-((self.h / self.scale) * row[j + 1]) * c0[source]);
                    sum.add(((self.h / self.scale) * row[j + 2]) * s1[source]);
                    sum.add(-((self.h / self.scale) * row[j + 3]) * c1[source]);
                    self.h * sum.value()
                };
                let (dv, df) = (clock(velocity), clock(force));
                work[source].add((self.work_scale * dv) * wf);
                work[source].add(wv * (self.work_scale * df));
                let term = (self.heat_scale * dv) * hv;
                heat[source].add(term);
                heat[source].add(term);
            }
        }
        let mut result = [ModalTonePowerDerivative::default(); MAX_HARMONIC_TONES];
        for source in 0..source_to_unique.len() {
            poll(cancelled)?;
            result[source] = finite(ModalTonePowerDerivative {
                input_work: work[source].value(),
                viscous_dissipation: heat[source].value(),
            })?;
        }
        poll(cancelled)?;
        Ok(result)
    }
}
