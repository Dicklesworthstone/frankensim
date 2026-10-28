//! Cache the existing exact quadratic flow integral for a fixed operator/clock.
//! Arbitrary signed quadratic observations include internal port work; none are
//! relabeled as dissipation or used to change the physical generator.
use super::*;
use crate::PreparedStepError;

fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl LinearSplit {
    /// Original declared storage matrix in the owner's unchanged state order.
    pub fn quadratic_matrix(&self) -> &[f64] {
        &self.quadratic.matrix
    }
}

/// Complete endpoint and time integrals from one held-input trajectory.
pub struct PreparedQuadraticStep {
    pub state: Vec<f64>,
    /// Quadratic forms first, followed by linear observations, in supplied order.
    pub integrals: Vec<f64>,
}

/// Fixed original state/port layout; no physical state is stored or reset here.
/// Warm work is matrix-vector multiplication and quadratic contractions only.
pub struct PreparedQuadraticFlow {
    n: usize,
    m: usize,
    d: usize,
    forward: Vec<f64>,
    inverse: Vec<f64>,
    transition: Vec<f64>,
    integrals: Vec<Vec<f64>>,
}
impl PreparedQuadraticFlow {
    /// Prepare the SAME energy-balanced exact-ZOH generator and stable Gramian
    /// doubling used by `step_render`. Inputs remain constant during a step.
    /// `quadratic_forms` contain original-state C for integral x^T C x dt;
    /// `linear_forms` contain original-state c for integral c^T x dt. Their
    /// values may be signed; the caller owns each observation's physical meaning.
    ///
    /// The explicit dimension covers states + held inputs + one constant
    /// coordinate and is capped at 256. Byte/work limits cover all caches and
    /// conservative 32-term/48-doubling preparation, before matrix allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        quadratic_forms: &[Vec<f64>],
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
    ) -> Result<Self, PhsError> {
        Self::new_controlled(
            sys,
            split,
            dt,
            quadratic_forms,
            linear_forms,
            budget,
            || false,
        )
        .map_err(|error| match error {
            PreparedStepError::Solver(error) => error,
            PreparedStepError::Cancelled => refuse("uncancelled quadratic preparation interrupted"),
        })
    }

    /// Same preparation with cooperative cancellation between transform rows,
    /// forms and after each existing Gramian call. One original bounded
    /// `flow_integral` is the largest atomic operation (at most 256 coordinates,
    /// 32 series terms and 48 doublings, within the supplied work/byte budgets).
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        quadratic_forms: &[Vec<f64>],
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        check_cancel(&mut cancelled)?;
        let (n, m) = (split.n, split.m);
        let d = n
            .checked_add(m)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| refuse("prepared quadratic dimension overflow"))?;
        let count = quadratic_forms
            .len()
            .checked_add(linear_forms.len())
            .ok_or_else(|| refuse("prepared observation count overflow"))?;
        let bytes = d
            .checked_mul(d)
            .and_then(|v| v.checked_mul(8))
            .and_then(|v| v.checked_mul(count.checked_mul(2)?.checked_add(24)?))
            .ok_or_else(|| refuse("prepared quadratic byte overflow"))?;
        let work = (d as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(224))
            .and_then(|v| v.checked_mul(count.max(1) as u64))
            .ok_or_else(|| refuse("prepared quadratic work overflow"))?;
        if d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || d > budget.max_dimension
            || bytes > budget.max_dense_bytes
            || work > budget.max_multiply_terms
        {
            return Err(refuse(
                "prepared quadratic complete operator exceeds explicit dimension/byte/work budget",
            )
            .into());
        }
        if !split.quadratic.exact
            || sys.n != n
            || sys.m != m
            || sys.j != split.j
            || sys.r != split.r
            || sys.g != split.g
            || sys.storage.quadratic_part().as_ref() != Some(&split.quadratic)
            || !dt.is_finite()
            || dt <= 0.0
            || quadratic_forms
                .iter()
                .any(|c| c.len() != n * n || !finite(c))
            || linear_forms.iter().any(|c| c.len() != n || !finite(c))
        {
            return Err(refuse(
                "prepared quadratic requires matching exact storage/structure and finite observations/clock",
            ).into());
        }
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            let weight = if split.values[i] > 0.0 {
                det::sqrt(split.values[i])
            } else {
                1.0
            };
            for j in 0..n {
                forward[i * n + j] = weight * split.basis[j * n + i];
                inverse[j * n + i] = split.basis[j * n + i] / weight;
            }
        }
        let mut flow = vec![0.0; n * n];
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            for k in 0..n {
                for j in 0..n {
                    flow[i * n + j] +=
                        (sys.j[i * n + k] - sys.r[i * n + k]) * split.quadratic.matrix[k * n + j];
                }
            }
        }
        check_cancel(&mut cancelled)?;
        let balanced = product(&product(&forward, &flow, n), &inverse, n);
        check_cancel(&mut cancelled)?;
        let mut generator = vec![0.0; d * d];
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            generator[i * d..i * d + n].copy_from_slice(&balanced[i * n..(i + 1) * n]);
            for p in 0..m {
                for j in 0..n {
                    generator[i * d + n + p] += forward[i * n + j] * sys.g[j * m + p];
                }
            }
        }
        let mut transpose = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                transpose[i * n + j] = inverse[j * n + i];
            }
        }
        let mut transition = Vec::new();
        let mut integrals = Vec::with_capacity(count);
        for index in 0..count.max(1) {
            check_cancel(&mut cancelled)?;
            let mut form = vec![0.0; d * d];
            if index < quadratic_forms.len() {
                let c = product(
                    &product(&transpose, &quadratic_forms[index], n),
                    &inverse,
                    n,
                );
                for i in 0..n {
                    form[i * d..i * d + n].copy_from_slice(&c[i * n..(i + 1) * n]);
                }
            } else if index < count {
                let c = &linear_forms[index - quadratic_forms.len()];
                for i in 0..n {
                    let value = 0.5 * (0..n).map(|j| c[j] * inverse[j * n + i]).sum::<f64>();
                    form[i * d + d - 1] = value;
                    form[(d - 1) * d + i] = value;
                }
            }
            check_cancel(&mut cancelled)?;
            let prepared = gramian::flow_integral(&generator, &form, d, dt).map_err(refuse)?;
            check_cancel(&mut cancelled)?;
            if transition.is_empty() {
                transition = prepared.transition;
            }
            if index < count {
                integrals.push(prepared.loss_integral);
            }
        }
        check_cancel(&mut cancelled)?;
        Ok(Self {
            n,
            m,
            d,
            forward,
            inverse,
            transition,
            integrals,
        })
    }

    /// Evaluate one complete candidate; accepted physical state belongs to the
    /// caller and is untouched on every failure. No matrix exponential is run.
    pub fn evaluate(
        &self,
        state: &[f64],
        input: &[f64],
    ) -> Result<PreparedQuadraticStep, PhsError> {
        if state.len() != self.n || input.len() != self.m || !finite(state) || !finite(input) {
            return Err(refuse(
                "prepared quadratic state/input layout or finiteness",
            ));
        }
        let mut z = vec![0.0; self.d];
        for i in 0..self.n {
            z[i] = (0..self.n)
                .map(|j| self.forward[i * self.n + j] * state[j])
                .sum();
        }
        z[self.n..self.n + self.m].copy_from_slice(input);
        z[self.d - 1] = 1.0;
        let mut balanced = vec![0.0; self.n];
        for i in 0..self.n {
            balanced[i] = (0..self.d)
                .map(|j| self.transition[i * self.d + j] * z[j])
                .sum();
        }
        let mut endpoint = vec![0.0; self.n];
        for i in 0..self.n {
            endpoint[i] = (0..self.n)
                .map(|j| self.inverse[i * self.n + j] * balanced[j])
                .sum();
        }
        let integrals: Vec<_> = self.integrals.iter().map(|w| quadratic(w, &z)).collect();
        if !finite(&endpoint) || !finite(&integrals) {
            return Err(refuse("prepared quadratic candidate overflows"));
        }
        Ok(PreparedQuadraticStep {
            state: endpoint,
            integrals,
        })
    }
}
