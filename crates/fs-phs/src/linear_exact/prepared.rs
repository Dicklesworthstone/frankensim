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
    /// Whether the original storage is exactly quadratic, including any
    /// retained linear material-memory states.
    pub fn is_exact_quadratic(&self) -> bool { self.quadratic.exact }

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
    dt: f64,
    growth_factor: f64,
    balanced_input: Vec<f64>,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            sys,
            split,
            dt,
            quadratic_forms,
            linear_forms,
            budget,
            cancelled,
            None,
        )
    }

    /// Append original linear states with numerical coordinate scales, retaining
    /// the base owner's energy transform and the SAME exponential/Gramian kernel.
    /// `generator_per_s` is the complete original-state generator, including the
    /// attachment reaction; inputs act through the unchanged base source ports.
    /// Positive `extra_scales` condition the added coordinates only: they confer
    /// no storage, passivity or dissipation meaning on an attached realization.
    /// The caller must retain its external work and original admission reports.
    #[allow(clippy::too_many_arguments)]
    pub fn new_extended_controlled(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        generator_per_s: &[f64],
        extra_scales: &[f64],
        dt: f64,
        quadratic_forms: &[Vec<f64>],
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            sys,
            split,
            dt,
            quadratic_forms,
            linear_forms,
            budget,
            cancelled,
            Some((generator_per_s, extra_scales)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        quadratic_forms: &[Vec<f64>],
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
        extension: Option<(&[f64], &[f64])>,
    ) -> Result<Self, PreparedStepError> {
        check_cancel(&mut cancelled)?;
        let (base_n, m) = (split.n, split.m);
        let n = base_n
            .checked_add(extension.map_or(0, |(_, scales)| scales.len()))
            .ok_or_else(|| refuse("prepared extended dimension overflow"))?;
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
            .and_then(|v| {
                v.checked_mul(count.checked_mul(2)?.checked_add(if extension.is_some() {
                    36
                } else {
                    24
                })?)
            })
            .ok_or_else(|| refuse("prepared quadratic byte overflow"))?;
        let work = (d as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(224))
            .and_then(|v| v.checked_mul(count.max(1) as u64))
            .and_then(|v| {
                v.checked_add(
                    (d as u64)
                        .checked_pow(3)?
                        .checked_mul(if extension.is_some() { 400 } else { 0 })?,
                )
            })
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
            || sys.n != base_n
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
            || extension.is_some_and(|(a, scales)| {
                a.len() != n * n || !finite(a) || scales.iter().any(|v| !v.is_finite() || *v <= 0.0)
            })
        {
            return Err(refuse(
                "prepared quadratic requires matching exact storage/structure and finite observations/clock",
            ).into());
        }
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..base_n {
            check_cancel(&mut cancelled)?;
            let weight = if split.values[i] > 0.0 {
                det::sqrt(split.values[i])
            } else {
                1.0
            };
            for j in 0..base_n {
                forward[i * n + j] = weight * split.basis[j * base_n + i];
                inverse[j * n + i] = split.basis[j * base_n + i] / weight;
            }
        }
        let flow = if let Some((generator, scales)) = extension {
            for (i, scale) in scales.iter().enumerate() {
                let index = base_n + i;
                forward[index * n + index] = *scale;
                inverse[index * n + index] = 1.0 / scale;
            }
            generator.to_vec()
        } else {
            let mut flow = vec![0.0; n * n];
            for i in 0..n {
                check_cancel(&mut cancelled)?;
                for k in 0..n {
                    for j in 0..n {
                        flow[i * n + j] += (sys.j[i * n + k] - sys.r[i * n + k])
                            * split.quadratic.matrix[k * n + j];
                    }
                }
            }
            flow
        };
        check_cancel(&mut cancelled)?;
        let balanced = product(&product(&forward, &flow, n), &inverse, n);
        check_cancel(&mut cancelled)?;
        // An extended external realization needs a sharp bound: Gershgorin
        // alone can interpret cancellation in a negative rank-one port damper
        // as enormous positive growth. Reuse the original bounded Jacobi owner.
        let mut growth_rate = 0.0_f64;
        if extension.is_some() {
            growth_rate = extended_growth_rate(&balanced, n, &mut cancelled)?;
        } else {
            for i in 0..n {
                let row = balanced[i * n + i]
                    + (0..n)
                        .filter(|j| *j != i)
                        .map(|j| (0.5 * balanced[i * n + j] + 0.5 * balanced[j * n + i]).abs())
                        .sum::<f64>();
                growth_rate = growth_rate.max(row);
            }
        }
        let growth_factor = det::exp(growth_rate * dt);
        // The original constructor does not need this envelope: a semidefinite
        // passive system may have finite polynomial flow but an overflowing
        // Euclidean exponential bound. Keep its prior admission unchanged.
        if !finite(&balanced) || (extension.is_some() && !growth_factor.is_finite()) {
            return Err(refuse("extended coordinate growth bound is not finite").into());
        }
        let mut generator = vec![0.0; d * d];
        let mut balanced_input = vec![0.0; n * m];
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            generator[i * d..i * d + n].copy_from_slice(&balanced[i * n..(i + 1) * n]);
            for p in 0..m {
                for j in 0..base_n {
                    generator[i * d + n + p] += forward[i * n + j] * sys.g[j * m + p];
                }
                balanced_input[i * m + p] = generator[i * d + n + p];
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
            dt,
            growth_factor,
            balanced_input,
        })
    }

    /// Bound the numerical coordinate norm throughout the complete interval:
    /// exp(max(0,mu)*dt) * (||T x0|| + dt*||T G u||). The retained base storage
    /// is at most half this norm squared, even when added states have NO known
    /// physical storage. No terminal energy or all-frequency passivity follows.
    pub fn interval_coordinate_norm_bound(
        &self,
        state: &[f64],
        input: &[f64],
    ) -> Result<f64, PhsError> {
        if state.len() != self.n || input.len() != self.m || !finite(state) || !finite(input) {
            return Err(refuse(
                "coordinate interval bound state/input layout or finiteness",
            ));
        }
        let mut initial = 0.0;
        let mut drive = 0.0;
        for i in 0..self.n {
            let x: f64 = (0..self.n)
                .map(|j| self.forward[i * self.n + j] * state[j])
                .sum();
            let f: f64 = (0..self.m)
                .map(|j| self.balanced_input[i * self.m + j] * input[j])
                .sum();
            initial += x * x;
            drive += f * f;
        }
        let bound = self.growth_factor * (det::sqrt(initial) + self.dt * det::sqrt(drive));
        if !bound.is_finite() {
            return Err(refuse("coordinate interval bound overflows"));
        }
        Ok(bound)
    }

    /// Evaluate one complete candidate; accepted physical state belongs to the
    /// caller and is untouched on every failure. No matrix exponential is run.
    pub fn evaluate(
        &self,
        state: &[f64],
        input: &[f64],
    ) -> Result<PreparedQuadraticStep, PhsError> {
        let mut z = vec![0.0; self.d];
        let mut balanced = vec![0.0; self.n];
        let mut endpoint = vec![0.0; self.n];
        let mut integrals = vec![0.0; self.integrals.len()];
        self.evaluate_into(state, input, &mut z, &mut balanced, &mut endpoint, &mut integrals)?;
        Ok(PreparedQuadraticStep { state: endpoint, integrals })
    }

    // Shared contractions for allocating evaluation and a prepared exact owner.
    // Output is numerical scratch; the public caller publishes only on success.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_into(
        &self, state: &[f64], input: &[f64], z: &mut [f64], balanced: &mut [f64],
        endpoint: &mut [f64], integrals: &mut [f64],
    ) -> Result<(), PhsError> {
        if state.len() != self.n || input.len() != self.m || !finite(state) || !finite(input)
            || z.len() != self.d || balanced.len() != self.n || endpoint.len() != self.n
            || integrals.len() != self.integrals.len()
        {
            return Err(refuse("prepared quadratic state/input/scratch layout or finiteness"));
        }
        for i in 0..self.n {
            z[i] = (0..self.n)
                .map(|j| self.forward[i * self.n + j] * state[j])
                .sum();
        }
        z[self.n..self.n + self.m].copy_from_slice(input);
        z[self.d - 1] = 1.0;
        for i in 0..self.n {
            balanced[i] = (0..self.d)
                .map(|j| self.transition[i * self.d + j] * z[j])
                .sum();
        }
        for i in 0..self.n {
            endpoint[i] = (0..self.n)
                .map(|j| self.inverse[i * self.n + j] * balanced[j])
                .sum();
        }
        for (out, w) in integrals.iter_mut().zip(&self.integrals) { *out = quadratic(w, z); }
        if !finite(endpoint) || !finite(integrals) {
            return Err(refuse("prepared quadratic candidate overflows"));
        }
        Ok(())
    }

}

// A floating-point envelope of lambda_max(sym(A)), with recomputed reconstruction
// and orthogonality allowances. This is numerical, not a storage certificate.
fn extended_growth_rate(
    a: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, PreparedStepError> {
    if !finite(a) {
        return Err(refuse("nonfinite extended generator").into());
    }
    admit_jacobi_eigh(n).map_err(|_| refuse("extended interval Jacobi admission"))?;
    let mut symmetric = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            symmetric[i * n + j] = 0.5 * a[i * n + j] + 0.5 * a[j * n + i];
        }
    }
    check_cancel(cancelled)?;
    let (values, vectors) = checked_eigh(&symmetric, n)?;
    check_cancel(cancelled)?;
    let mut scaled = vectors.clone();
    let mut transpose = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            scaled[i * n + j] *= values[j];
            transpose[i * n + j] = vectors[j * n + i];
        }
    }
    let reconstructed = product(&scaled, &transpose, n);
    check_cancel(cancelled)?;
    let gram = product(&transpose, &vectors, n);
    let norm = |matrix: &[f64]| {
        matrix
            .chunks_exact(n)
            .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max)
    };
    let mut residual = 0.0_f64;
    let mut defect = 0.0_f64;
    for i in 0..n {
        check_cancel(cancelled)?;
        residual = residual.max(
            (0..n)
                .map(|j| (symmetric[i * n + j] - reconstructed[i * n + j]).abs())
                .sum(),
        );
        defect = defect.max(
            (0..n)
                .map(|j| (gram[i * n + j] - if i == j { 1.0 } else { 0.0 }).abs())
                .sum(),
        );
    }
    // Accumulation-length and magnitude allowance for both dense products,
    // subtraction and row norms; cancellation cannot erase their error scale.
    let roundoff = 64.0 * n as f64 * f64::EPSILON;
    let vector_product = norm(&vectors) * norm(&transpose);
    let lambda_magnitude = values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    residual += roundoff * (norm(&symmetric) + vector_product * lambda_magnitude);
    defect += roundoff * (1.0 + vector_product);
    let positive = values.iter().copied().fold(0.0_f64, f64::max);
    // S = V Lambda V^T + E. Its largest eigenvalue is bounded above by
    // max(0,max Lambda)*||V||_2^2 + ||E||_2, with ||V||_2^2 <= 1+defect.
    let bound = positive * (1.0 + defect) + residual;
    if !finite(&reconstructed) || !finite(&gram) || !bound.is_finite() {
        return Err(refuse("extended interval reconstruction bound is not finite").into());
    }
    Ok(bound)
}
