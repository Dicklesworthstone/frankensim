//! One-way observation of the complete continuous physical source trajectory.
use super::{PreparedExactTones, Sum, combine_tones, uncancelled};
use crate::linear_exact::{
    DenseExponentialBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, exponential_with_budget, finite,
    refuse,
};
use crate::{PhsError, PreparedStepError};
use std::sync::Arc;

/// Fixed-clock causal observation `z'=A*z+B*x`, `y=C*z+D*x` of the original
/// complete physical state x under continuous tone forcing. Matrices use the
/// caller's original coordinates, not energy-normalized/modal replacements.
///
/// Preparation reuses the source's exact balanced generator and the existing
/// bounded matrix exponential. Only observer endpoint and output rows are
/// cached. Observer states carry no physical storage, heat, loading or reaction;
/// the source's original accepted step remains the authority for all physics.
/// A finite realization is not certified passive, stable or bandlimited here.
///
/// The cache is bound to this exact source owner, including clock and source
/// ordering. Moving the owner is safe; an independently rebuilt owner refuses.
/// Physical and observer states remain caller-owned and checkpointable.
pub struct PreparedToneObserver {
    identity: Arc<()>,
    states: usize,
    outputs: usize,
    width: usize,
    state_rows: Vec<f64>,
    output_rows: Vec<f64>,
    coordinates: Vec<f64>,
    input: Vec<f64>,
    sums: Vec<Sum>,
    candidate_state: Vec<f64>,
    candidate_output: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl PreparedToneObserver {
    /// Prepare row-major A(K×K), B(K×N), C(O×K), D(O×N). N includes every
    /// original physical coordinate. K may be zero for pure direct feedthrough;
    /// at least one output is required. Output is evaluated at the interval end.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source: &PreparedExactTones<'_>,
        observer_dim: usize,
        output_dim: usize,
        a: &[f64],
        b: &[f64],
        c: &[f64],
        d: &[f64],
        budget: DenseExponentialBudget,
    ) -> Result<Self, PhsError> {
        Self::new_controlled(source, observer_dim, output_dim, a, b, c, d, budget, || {
            false
        })
        .map_err(uncancelled)
    }

    /// Same cold preparation with polling between rows and before/after the
    /// existing bounded exponential (at most 256 complete coordinates). No
    /// source state, accepted ledger or partially prepared cache is published.
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        source: &PreparedExactTones<'_>,
        observer_dim: usize,
        output_dim: usize,
        a: &[f64],
        b: &[f64],
        c: &[f64],
        d: &[f64],
        budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let (n, k, o) = (source.state_dim(), observer_dim, output_dim);
        let s = source.flow.source_scales().len();
        let base = n + s;
        let width = base
            .checked_add(k)
            .ok_or_else(|| refuse("tone observer dimension overflow"))?;
        if o == 0
            || width > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || width > budget.max_dimension
            || k.checked_mul(k) != Some(a.len())
            || k.checked_mul(n) != Some(b.len())
            || o.checked_mul(k) != Some(c.len())
            || o.checked_mul(n) != Some(d.len())
            || [a, b, c, d].iter().any(|matrix| !finite(matrix))
        {
            return Err(refuse(
                "tone observer requires complete finite matrices within dimension budget",
            )
            .into());
        }
        // Includes the complete augmented exponential workspace, retained rows,
        // candidates and compensated source packing before any dense allocation.
        let (w, kr, or, sr) = (width as u128, k as u128, o as u128, s as u128);
        let bytes = 8 * (6 * w * w + (kr + or) * w + w + 2 * (kr + or) + 3 * sr);
        let work = 80 * w * w * w + or * (n as u128 + kr) * w + w * w;
        if bytes > budget.max_dense_bytes as u128 || work > budget.max_multiply_terms as u128 {
            return Err(refuse("complete tone observer exceeds explicit byte/work budget").into());
        }
        let (physical, transform) = source
            .flow
            .tone_operator()
            .ok_or_else(|| refuse("tone observer requires original continuous source operator"))?;
        let mut operator = vec![0.0; width * width];
        for i in 0..base {
            poll(&mut cancelled)?;
            for j in 0..base {
                operator[i * width + j] = physical[i * (base + 1) + j] * source.duration;
            }
        }
        for i in 0..k {
            poll(&mut cancelled)?;
            for j in 0..n {
                operator[(base + i) * width + j] =
                    (b[i * n + j] / transform[j * n + j]) * source.duration;
            }
            for j in 0..k {
                operator[(base + i) * width + base + j] = a[i * k + j] * source.duration;
            }
        }
        poll(&mut cancelled)?;
        let transition = exponential_with_budget(&operator, width, budget)?;
        poll(&mut cancelled)?;
        let mut state_rows = vec![0.0; k * width];
        for i in 0..k {
            state_rows[i * width..(i + 1) * width]
                .copy_from_slice(&transition[(base + i) * width..(base + i + 1) * width]);
        }
        let mut output_rows = vec![0.0; o * width];
        for i in 0..o {
            poll(&mut cancelled)?;
            for j in 0..width {
                let mut value = Sum::default();
                for p in 0..k {
                    value.add(c[i * k + p] * transition[(base + p) * width + j]);
                }
                for p in 0..n {
                    value.add((d[i * n + p] / transform[p * n + p]) * transition[p * width + j]);
                }
                output_rows[i * width + j] = value.value();
            }
        }
        if !finite(&output_rows) {
            return Err(refuse("tone observer output rows overflow").into());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            identity: Arc::clone(&source.identity),
            states: k,
            outputs: o,
            width,
            state_rows,
            output_rows,
            coordinates: vec![0.0; width],
            input: vec![0.0; s],
            sums: vec![Sum::default(); s],
            candidate_state: vec![0.0; k],
            candidate_output: vec![0.0; o],
        })
    }

    /// Number of original caller-defined observation states.
    pub fn observer_state_dim(&self) -> usize {
        self.states
    }
    /// Number of endpoint observation outputs.
    pub fn output_dim(&self) -> usize {
        self.outputs
    }

    /// Preview the exact observer endpoint and output on the source's complete
    /// within-interval trajectory. Source coefficients have the same port-major
    /// caller ordering as its physical step, including exact duplicate tones.
    /// No allocation, exponential, transcendental or physical step occurs here.
    /// Both caller outputs remain unchanged on every refusal. The consumer must
    /// publish them only alongside an accepted physical step using these inputs.
    #[allow(clippy::too_many_arguments)]
    pub fn preview_into(
        &mut self,
        source: &PreparedExactTones<'_>,
        physical_initial: &[f64],
        observer_initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        observer_next: &mut [f64],
        output: &mut [f64],
    ) -> Result<(), PhsError> {
        self.preview_into_controlled(
            source,
            physical_initial,
            observer_initial,
            cosine,
            sine,
            observer_next,
            output,
            || false,
        )
        .map_err(uncancelled)
    }

    /// Same allocation-free candidate preview with cooperative cancellation.
    /// Scratch may change; caller state, output and the bound source do not.
    #[allow(clippy::too_many_arguments)]
    pub fn preview_into_controlled(
        &mut self,
        source: &PreparedExactTones<'_>,
        physical_initial: &[f64],
        observer_initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        observer_next: &mut [f64],
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), PreparedStepError> {
        self.preview_impl(
            source,
            physical_initial,
            observer_initial,
            cosine,
            sine,
            None,
            observer_next,
            output,
            cancelled,
        )
    }

    /// Observe affine tone envelopes through the SAME continuous source and
    /// filter map. Four endpoint arrays use the interval-start carrier phase,
    /// exactly as `PreparedExactTones::step_affine_into`. Requires affine source
    /// preparation; physical and observation candidates still commit together.
    #[allow(clippy::too_many_arguments)]
    pub fn preview_affine_into(
        &mut self,
        source: &PreparedExactTones<'_>,
        physical_initial: &[f64],
        observer_initial: &[f64],
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        observer_next: &mut [f64],
        output: &mut [f64],
    ) -> Result<(), PhsError> {
        self.preview_affine_into_controlled(
            source,
            physical_initial,
            observer_initial,
            cosine_start,
            sine_start,
            cosine_end,
            sine_end,
            observer_next,
            output,
            || false,
        )
        .map_err(uncancelled)
    }

    /// Same allocation-free affine candidate preview with cancellation. All
    /// caller outputs remain unchanged on any failure or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn preview_affine_into_controlled(
        &mut self,
        source: &PreparedExactTones<'_>,
        physical_initial: &[f64],
        observer_initial: &[f64],
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        observer_next: &mut [f64],
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), PreparedStepError> {
        self.preview_impl(
            source,
            physical_initial,
            observer_initial,
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
            observer_next,
            output,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn preview_impl(
        &mut self,
        source: &PreparedExactTones<'_>,
        physical_initial: &[f64],
        observer_initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        end: Option<(&[f64], &[f64])>,
        observer_next: &mut [f64],
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), PreparedStepError> {
        poll(&mut cancelled)?;
        let n = source.state_dim();
        let s = self.input.len();
        if !Arc::ptr_eq(&self.identity, &source.identity)
            || source.is_affine() != end.is_some()
            || physical_initial.len() != n
            || observer_initial.len() != self.states
            || cosine.len() != source.port_dim() * source.tone_count()
            || sine.len() != cosine.len()
            || observer_next.len() != self.states
            || output.len() != self.outputs
            || end.is_some_and(|(c, s)| {
                c.len() != cosine.len() || s.len() != sine.len() || !finite(c) || !finite(s)
            })
            || [physical_initial, observer_initial, cosine, sine]
                .iter()
                .any(|v| !finite(v))
        {
            return Err(refuse(
                "tone observer source binding, state/input/output layout or finiteness",
            )
            .into());
        }
        if let Some((end_cosine, end_sine)) = end {
            super::affine::combine(
                &source.mapping,
                source.unique.len(),
                cosine,
                sine,
                end_cosine,
                end_sine,
                &mut self.sums,
                &mut self.input,
            )?;
        } else {
            combine_tones(
                &source.mapping,
                source.unique.len(),
                cosine,
                sine,
                &mut self.sums,
                &mut self.input,
            )?;
        }
        let (_, transform) = source
            .flow
            .tone_operator()
            .ok_or_else(|| refuse("tone observer source operator missing"))?;
        for i in 0..n {
            self.coordinates[i] = transform[i * n + i] * physical_initial[i];
        }
        for i in 0..s {
            self.coordinates[n + i] = source.flow.source_scales()[i] * self.input[i];
        }
        self.coordinates[n + s..].copy_from_slice(observer_initial);
        if !finite(&self.coordinates) {
            return Err(refuse("tone observer initial coordinates overflow").into());
        }
        for (rows, candidate) in [
            (&self.state_rows, &mut self.candidate_state),
            (&self.output_rows, &mut self.candidate_output),
        ] {
            for (row, out) in rows.chunks_exact(self.width).zip(candidate) {
                poll(&mut cancelled)?;
                let mut value = Sum::default();
                for (&coefficient, &coordinate) in row.iter().zip(&self.coordinates) {
                    value.add(coefficient * coordinate);
                }
                *out = value.value();
                if !out.is_finite() {
                    return Err(refuse("tone observer endpoint/output overflow").into());
                }
            }
        }
        poll(&mut cancelled)?;
        observer_next.copy_from_slice(&self.candidate_state);
        output.copy_from_slice(&self.candidate_output);
        Ok(())
    }
}
