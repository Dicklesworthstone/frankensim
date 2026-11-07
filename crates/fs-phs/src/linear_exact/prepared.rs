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

enum CachedIntegral {
    // Full quadratic observations retain their original property check and
    // contraction, including matrices that happen to be exactly affine.
    Dense { matrix: Vec<f64>, affine: bool },
    // Symmetric final row/column of a strictly linear observation's Gramian.
    Linear(Vec<f64>),
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
    integrals: Vec<CachedIntegral>,
    dt: f64,
    growth_factor: f64,
    storage_norm_factor: f64,
    balanced_input: Vec<f64>,
    source_scales: Vec<f64>,
    rotating_sources: bool,
    tone_generator: Option<Vec<f64>>,
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
    /// forms, linear Taylor stages and convergence-group doublings. One original
    /// full `flow_integral` remains the largest atomic operation (at most 256
    /// coordinates, 32 series terms and 48 doublings, within the supplied budgets).
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
            None,
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
            None,
            None,
        )
    }

    // The same complete operator/Gramian owner, with continuous source
    // quadratures instead of held source coordinates. All physical states and
    // original linear observations remain in their existing order.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_tones_controlled(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        frequencies: &[f64],
        affine: bool,
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            sys, split, dt, &[], linear_forms, budget, cancelled, None,
            Some((frequencies, affine, sys.m)), None,
        )
    }

    // Original physical ports rotate; the remaining shadow ports retain one
    // held coordinate each. Nonlinear discrete efforts are never replicated
    // across frequencies or assigned fictitious source/storage energy.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_mixed_tones_controlled(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        frequencies: &[f64],
        affine: bool,
        rotating_ports: usize,
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            sys, split, dt, &[], linear_forms, budget, cancelled, None,
            Some((frequencies, affine, rotating_ports)), None,
        )
    }

    pub(super) fn source_scales(&self) -> &[f64] {
        &self.source_scales
    }

    // Exact prepared continuous generator and original diagonal state map.
    // One-way observers must reuse these, including every source quadrature.
    pub(super) fn tone_operator(&self) -> Option<(&[f64], &[f64])> {
        self.tone_generator.as_deref().map(|generator| (generator, self.forward.as_slice()))
    }

    /// Retain the same extended generator and physical base storage while
    /// supplying the complete original-coordinate input map. This supports
    /// force-record states whose derivatives, rather than the mechanical
    /// coordinates, receive the held input. The input count remains the base
    /// system's port count. Input units and work forms belong to the caller;
    /// no passivity or physical meaning is inferred for the added coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn new_extended_input_controlled(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        generator_per_s: &[f64],
        extra_scales: &[f64],
        input_map: &[f64],
        dt: f64,
        quadratic_forms: &[Vec<f64>],
        linear_forms: &[Vec<f64>],
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(sys, split, dt, quadratic_forms, linear_forms, budget,
            cancelled, Some((generator_per_s, extra_scales)), None, Some(input_map))
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
        tones: Option<(&[f64], bool, usize)>,
        input_override: Option<&[f64]>,
    ) -> Result<Self, PreparedStepError> {
        check_cancel(&mut cancelled)?;
        let (base_n, ports) = (split.n, split.m);
        let m = if let Some((frequencies, affine, rotating_ports)) = tones {
            if frequencies.is_empty()
                || rotating_ports == 0
                || rotating_ports > ports
                || extension.is_some()
                || input_override.is_some()
                || frequencies.iter().any(|v| !v.is_finite() || *v < 0.0)
                || (affine && !(1.0 / dt).is_finite())
            {
                return Err(refuse(
                    "tone sources require nonempty finite frequencies and physical ports",
                ).into());
            }
            rotating_ports.checked_mul(frequencies.len())
                .and_then(|v| v.checked_mul(if affine { 4 } else { 2 }))
                .and_then(|v| v.checked_add(ports - rotating_ports))
                .ok_or_else(|| refuse("tone source dimension overflow"))?
        } else {
            ports
        };
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
                } else if tones.is_some() {
                    25
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
            || sys.m != ports
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
            || input_override.is_some_and(|input| input.len() != n * m || !finite(input))
        {
            return Err(refuse(
                "prepared quadratic requires matching exact storage/structure and finite observations/clock",
            ).into());
        }
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..base_n {
            check_cancel(&mut cancelled)?;
            if extension.is_some() {
                // External-state trajectory admission uses the original
                // energy-coordinate logarithmic growth envelope. Preserve
                // that metric: diagonal coordinates can give a very loose
                // positive growth bound for a stiff passive correlated Q.
                let scale = if split.values[i] > 0.0 {
                    det::sqrt(split.values[i])
                } else {
                    1.0
                };
                for j in 0..base_n {
                    forward[i * n + j] = scale * split.basis[j * base_n + i];
                    inverse[j * n + i] = split.basis[j * base_n + i] / scale;
                }
            } else {
                let scale = split.coordinate_scale(i);
                forward[i * n + i] = scale;
                inverse[i * n + i] = 1.0 / scale;
            }
        }
        // Diagonal equilibration does not whiten correlated storage. The
        // symmetric row-sum bound controls lambda_max(D^-T Q D^-1), so the
        // interval norm still bounds the ORIGINAL physical base energy.
        // Keeping at least one also bounds numerical extension coordinates.
        let mut storage_norm = 1.0_f64;
        let spectral_effort = if extension.is_some() {
            // These maps use the actual retained inverse, including its
            // expanded stride. Absolute products bound cancellation roundoff.
            let mut effort = vec![0.0; base_n * base_n];
            let mut absolute = vec![0.0; base_n * base_n];
            for i in 0..base_n {
                check_cancel(&mut cancelled)?;
                for j in 0..base_n {
                    for k in 0..base_n {
                        let term = split.quadratic.matrix[i * base_n + k]
                            * inverse[k * n + j];
                        effort[i * base_n + j] += term;
                        absolute[i * base_n + j] += term.abs();
                    }
                }
            }
            Some((effort, absolute))
        } else {
            None
        };
        for i in 0..base_n {
            check_cancel(&mut cancelled)?;
            let row: f64 = (0..base_n)
                .map(|j| {
                    if let Some((effort, absolute)) = &spectral_effort {
                        let mut value = 0.0;
                        let mut products = 0.0;
                        for k in 0..base_n {
                            value += inverse[k * n + i] * effort[k * base_n + j];
                            products += inverse[k * n + i].abs() * absolute[k * base_n + j];
                        }
                        value.abs() + 16.0 * base_n as f64 * f64::EPSILON * products
                    } else {
                        (split.quadratic.matrix[i * base_n + j]
                            / forward[i * n + i]
                            / forward[j * n + j])
                            .abs()
                    }
                })
                .sum();
            if !row.is_finite() {
                return Err(refuse("prepared storage norm row is not finite").into());
            }
            storage_norm = storage_norm.max(row);
        }
        // Numerical roundoff allowance for the two divisions and row sum;
        // this changes only the conservative envelope, never Q or the flow.
        let storage_norm_factor =
            det::sqrt(storage_norm * (1.0 + 8.0 * base_n as f64 * f64::EPSILON));
        if !storage_norm_factor.is_finite() {
            return Err(refuse("prepared storage norm bound is not finite").into());
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
        let mut source_scales = vec![1.0; m];
        let mut physical_input = vec![0.0; n * ports];
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            generator[i * d..i * d + n].copy_from_slice(&balanced[i * n..(i + 1) * n]);
            for p in 0..ports {
                if let Some(input) = input_override {
                    for j in 0..n {
                        physical_input[i * ports + p] += forward[i * n + j] * input[j * ports + p];
                    }
                } else {
                    for j in 0..base_n {
                        physical_input[i * ports + p] += forward[i * n + j] * sys.g[j * ports + p];
                    }
                }
            }
        }
        if let Some((frequencies, affine, rotating_ports)) = tones {
            let stride = if affine { 4 } else { 2 };
            for p in 0..rotating_ports {
                let magnitude = (0..n)
                    .map(|i| physical_input[i * ports + p].abs())
                    .fold(0.0_f64, f64::max);
                let scale = if magnitude == 0.0 { 1.0 } else { magnitude };
                if !scale.is_finite() {
                    return Err(refuse("tone input scale must be finite").into());
                }
                for (tone, &omega) in frequencies.iter().enumerate() {
                    let c = stride * (p * frequencies.len() + tone);
                    source_scales[c..c + stride].fill(scale);
                    generator[(n + c) * d + n + c + 1] = omega;
                    generator[(n + c + 1) * d + n + c] = -omega;
                    if affine {
                        // Delta coefficients have force units. Their rotating
                        // pair drives an affine envelope over this exact h.
                        generator[(n + c) * d + n + c + 2] = 1.0 / dt;
                        generator[(n + c + 1) * d + n + c + 3] = 1.0 / dt;
                        generator[(n + c + 2) * d + n + c + 3] = omega;
                        generator[(n + c + 3) * d + n + c + 2] = -omega;
                    }
                    for i in 0..n {
                        generator[i * d + n + c] = physical_input[i * ports + p] / scale;
                        balanced_input[i * m + c] = generator[i * d + n + c];
                    }
                }
            }
            let rotating = stride * rotating_ports * frequencies.len();
            for p in rotating_ports..ports {
                let c = rotating + p - rotating_ports;
                for i in 0..n {
                    generator[i * d + n + c] = physical_input[i * ports + p];
                    balanced_input[i * m + c] = physical_input[i * ports + p];
                }
            }
        } else {
            for i in 0..n {
                for p in 0..m {
                    generator[i * d + n + p] = physical_input[i * m + p];
                    balanced_input[i * m + p] = physical_input[i * m + p];
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
        // Quadratic observations retain the unchanged full Gramian path.
        // With no observations, one zero form still prepares the transition.
        let full_count = if count == 0 { 1 } else { quadratic_forms.len() };
        for index in 0..full_count {
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
            }
            check_cancel(&mut cancelled)?;
            let prepared = gramian::flow_integral(&generator, &form, d, dt).map_err(refuse)?;
            check_cancel(&mut cancelled)?;
            if transition.is_empty() {
                transition = prepared.transition;
            }
            if index < quadratic_forms.len() {
                let affine = is_affine_integral(&prepared.loss_integral, d);
                integrals.push(CachedIntegral::Dense {
                    matrix: prepared.loss_integral,
                    affine,
                });
            }
        }
        if !linear_forms.is_empty() {
            let mut rows = Vec::with_capacity(linear_forms.len());
            for c in linear_forms {
                check_cancel(&mut cancelled)?;
                let mut row = vec![0.0; d];
                for i in 0..n {
                    row[i] = 0.5 * (0..n).map(|j| c[j] * inverse[j * n + i]).sum::<f64>();
                }
                rows.push(row);
            }
            // The existing conservative dense-reference budgets dominate this
            // shared transition, one live group snapshot and compact row work.
            let prepared = gramian::affine_flow_integrals(
                &generator, rows, d, dt, &mut cancelled,
            ).map_err(refuse)?.ok_or(PreparedStepError::Cancelled)?;
            if transition.is_empty() {
                transition = prepared.transition;
            }
            integrals.extend(prepared.half_rows.into_iter().map(CachedIntegral::Linear));
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
            storage_norm_factor,
            balanced_input,
            source_scales,
            rotating_sources: tones.is_some(),
            tone_generator: tones.map(|_| generator),
        })
    }

    /// Bound the numerical coordinate norm throughout the complete interval:
    /// sqrt(max(1,rho)) * exp(max(0,mu)*dt) * (||T x0|| + dt*||T G u||),
    /// where rho bounds the original storage in the equilibrated coordinates.
    /// The retained base storage is at most half this norm squared, even when
    /// added states have NO known physical storage. No terminal energy or
    /// all-frequency passivity follows.
    pub fn interval_coordinate_norm_bound(
        &self,
        state: &[f64],
        input: &[f64],
    ) -> Result<f64, PhsError> {
        if self.rotating_sources {
            return Err(refuse(
                "held-input interval envelope is not a continuous-tone envelope",
            ));
        }
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
        let bound = self.storage_norm_factor
            * self.growth_factor
            * (det::sqrt(initial) + self.dt * det::sqrt(drive));
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

    /// Recover a declared linear observation's cached integral in ORIGINAL
    /// coordinates as `(state_coefficients, held_input_coefficients)`.
    /// The index follows `evaluate().integrals`, including any preceding
    /// quadratic forms. A nonlinear observation is refused. Its Gramian has
    /// only the constant row/column: no trial physical state, constitutive
    /// probe or second exponential is needed. This allocates a cold-path map;
    /// consumers may retain it with this immutable prepared owner.
    pub fn linear_integral_map(
        &self,
        index: usize,
    ) -> Result<(Vec<f64>, Vec<f64>), PhsError> {
        let form = self.integrals.get(index)
            .ok_or_else(|| refuse("prepared linear observation index"))?;
        let last = self.d - 1;
        let linear = match form {
            CachedIntegral::Linear(row) => row[last] == 0.0,
            CachedIntegral::Dense { matrix, .. } => {
                matrix[last * self.d + last] == 0.0
                    && (0..last).all(|i| (0..last).all(|j| matrix[i * self.d + j] == 0.0))
            }
        };
        if !linear {
            return Err(refuse("prepared observation is not a linear form"));
        }
        let coefficient = |i: usize| match form {
            CachedIntegral::Linear(row) => row[i] + row[i],
            CachedIntegral::Dense { matrix, .. } => {
                matrix[i * self.d + last] + matrix[last * self.d + i]
            }
        };
        let balanced: Vec<f64> = (0..self.n).map(coefficient).collect();
        let state: Vec<f64> = (0..self.n)
            .map(|j| (0..self.n).map(|i| balanced[i] * self.forward[i * self.n + j]).sum())
            .collect();
        let input: Vec<f64> = (0..self.m)
            .map(|p| coefficient(self.n + p) * self.source_scales[p])
            .collect();
        if !finite(&state) || !finite(&input) {
            return Err(refuse("prepared linear observation map is nonfinite"));
        }
        Ok((state, input))
    }

    /// Shared allocation-free contractions using caller-owned scratch.
    /// `z` has states + inputs + one entries; `balanced` and `endpoint` cover
    /// every state; `integrals` covers the declared quadratic and linear forms.
    /// Outputs are candidate scratch and can change on refusal. The caller
    /// admits physical budgets and publishes accepted state only after success.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_into(
        &self, state: &[f64], input: &[f64], z: &mut [f64], balanced: &mut [f64],
        endpoint: &mut [f64], integrals: &mut [f64],
    ) -> Result<(), PhsError> {
        self.evaluate_impl(state, input, z, balanced, endpoint, integrals, false,
            |_, _| Ok(())).map(|_| ())
    }

    /// Held physical-port observations retain their integrated low words
    /// through multiplication by the conjugate inputs. General observations
    /// and rotating sources keep the public evaluator's existing semantics.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_port_work_into(
        &self, state: &[f64], input: &[f64], z: &mut [f64], balanced: &mut [f64],
        endpoint: &mut [f64], integrals: &mut [f64],
    ) -> Result<f64, PhsError> {
        self.evaluate_port_work_observed_into(state, input, [], z, balanced, endpoint, integrals)
            .map(|(work, _)| work)
    }

    /// Contract every work contribution with the integrated port low words
    /// before rounding the public observations. All contributions see the same
    /// physical trajectory; none is recovered by subtracting two total works.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_port_work_observed_into<const K: usize>(
        &self, state: &[f64], input: &[f64], work_inputs: [&[f64]; K],
        z: &mut [f64], balanced: &mut [f64], endpoint: &mut [f64], integrals: &mut [f64],
    ) -> Result<(f64, [f64; K]), PhsError> {
        if self.rotating_sources || self.integrals.len() != self.m
            || self.integrals.iter().any(|form| !matches!(form, CachedIntegral::Linear(_)))
        {
            return Err(refuse("held port work requires one linear observation per input"));
        }
        if work_inputs.iter().any(|values| values.len() != self.m || !finite(values)) {
            return Err(refuse("held work observations require complete finite inputs"));
        }
        let mut observed = [Sum::default(); K];
        let supplied = self.evaluate_impl(state, input, z, balanced, endpoint, integrals, true,
            |port, integral| {
                for (work, values) in observed.iter_mut().zip(work_inputs) {
                    work.add_scaled(integral, values[port])
                        .ok_or_else(|| refuse("prepared observed work is unrepresentable"))?;
                }
                Ok(())
            })?;
        let mut work = [0.0; K];
        for (value, sum) in work.iter_mut().zip(observed) {
            *value = sum.value()
                .ok_or_else(|| refuse("prepared observed work is unrepresentable"))?;
        }
        Ok((supplied, work))
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_impl(
        &self, state: &[f64], input: &[f64], z: &mut [f64], balanced: &mut [f64],
        endpoint: &mut [f64], integrals: &mut [f64], physical_work: bool,
        mut observe: impl FnMut(usize, Sum) -> Result<(), PhsError>,
    ) -> Result<f64, PhsError> {
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
        for p in 0..self.m {
            z[self.n + p] = self.source_scales[p] * input[p];
        }
        z[self.d - 1] = 1.0;
        for i in 0..self.n {
            balanced[i] = crate::port_input_row(
                &self.transition[i * self.d..(i + 1) * self.d], |j| z[j],
            ).and_then(Sum::value)
                .ok_or_else(|| refuse("prepared transition contraction is unrepresentable"))?;
        }
        for i in 0..self.n {
            endpoint[i] = (0..self.n)
                .map(|j| self.inverse[i * self.n + j] * balanced[j])
                .sum();
        }
        let mut work = Sum::default();
        for (p, (out, form)) in integrals.iter_mut().zip(&self.integrals).enumerate() {
            *out = match form {
                CachedIntegral::Dense { matrix, affine: true } => affine_integral(matrix, z),
                CachedIntegral::Dense { matrix, affine: false } => quadratic(matrix, z),
                CachedIntegral::Linear(row) if physical_work => {
                    // z[last]=1: the off-diagonal half-row occurs twice;
                    // its final diagonal contributes only once.
                    let last = self.d - 1;
                    let half = crate::port_input_row(&row[..last], |j| z[j])
                        .ok_or_else(|| refuse("prepared port integral is unrepresentable"))?;
                    let mut integral = half;
                    integral.add_scaled(half, 1.0)
                        .and_then(|()| integral.product(row[last], 1.0))
                        .ok_or_else(|| refuse("prepared port integral is unrepresentable"))?;
                    work.add_scaled(integral, input[p])
                        .ok_or_else(|| refuse("prepared supplied work is unrepresentable"))?;
                    observe(p, integral)?;
                    integral.value()
                        .ok_or_else(|| refuse("prepared port integral is unrepresentable"))?
                }
                CachedIntegral::Linear(row) => compact_affine_integral(row, z),
            };
        }
        if !finite(endpoint) || !finite(integrals) {
            return Err(refuse("prepared quadratic candidate overflows"));
        }
        work.value().ok_or_else(|| refuse("prepared supplied work is unrepresentable"))
    }

}

fn is_affine_integral(matrix: &[f64], dimension: usize) -> bool {
    let last = dimension - 1;
    (0..last).all(|i| {
        matrix[i * dimension..i * dimension + last]
            .iter()
            .all(|x| *x == 0.0)
    })
}

// Same left-associated products and row-major additions as quadratic(). With
// finite coordinates, skipped products are signed zero; an accumulator that
// starts at +0 cannot change by adding them. Do not combine the two symmetric
// entries or contract a mean-output map: either would change the rounding.
fn affine_integral(matrix: &[f64], coordinates: &[f64]) -> f64 {
    let dimension = coordinates.len();
    let last = dimension - 1;
    let mut result = 0.0;
    for i in 0..last {
        result += coordinates[i] * matrix[i * dimension + last] * coordinates[last];
    }
    for j in 0..dimension {
        result += coordinates[last] * matrix[last * dimension + j] * coordinates[j];
    }
    result
}

// The same two half-row contractions as affine_integral(), with the symmetric
// entries stored once. Combining the halves would change cancellation rounding.
fn compact_affine_integral(row: &[f64], coordinates: &[f64]) -> f64 {
    let last = coordinates.len() - 1;
    let mut result = 0.0;
    for i in 0..last {
        result += coordinates[i] * row[i] * coordinates[last];
    }
    for j in 0..coordinates.len() {
        result += coordinates[last] * row[j] * coordinates[j];
    }
    result
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

#[cfg(test)]
mod affine_integral_tests {
    use super::*;

    fn dense(form: &CachedIntegral, dimension: usize) -> Vec<f64> {
        match form {
            CachedIntegral::Dense { matrix, .. } => matrix.clone(),
            CachedIntegral::Linear(row) => {
                let mut matrix = vec![0.0; dimension * dimension];
                let last = dimension - 1;
                for i in 0..dimension {
                    matrix[i * dimension + last] = row[i];
                    matrix[last * dimension + i] = row[i];
                }
                matrix
            }
        }
    }

    #[test]
    fn affine_observations_preserve_dense_bits_without_discarding_quadratic_terms() {
        let system = crate::mass_spring_damper(1.0, 9.0, 0.7).unwrap();
        let split = super::super::split_quadratic(
            &system,
            super::super::LinearSplitBudget {
                max_modes: 2,
                max_omega_dt: 2.8,
            },
        )
        .unwrap();
        let flow = PreparedQuadraticFlow::new(
            &system,
            &split,
            0.03,
            &[vec![0.0, 0.0, 0.0, 1.0]],
            &[vec![0.3, -0.2]],
            DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX,
                max_multiply_terms: u64::MAX,
            },
        )
        .unwrap();
        assert!(matches!(flow.integrals[0], CachedIntegral::Dense { affine: false, .. }));
        assert!(matches!(flow.integrals[1], CachedIntegral::Linear(_)));
        for scale in [1.0 / 1_099_511_627_776.0, 1.0, 1_099_511_627_776.0] {
            for i in 0..32 {
                let state = [
                    scale * (i as f64 - 15.0) / 13.0,
                    scale * (7.0 - i as f64) / 19.0,
                ];
                let input = [scale * (i as f64 % 7.0 - 3.0)];
                let mut z = vec![0.0; flow.d];
                let mut values = [0.0; 2];
                flow.evaluate_into(
                    &state,
                    &input,
                    &mut z,
                    &mut [0.0; 2],
                    &mut [0.0; 2],
                    &mut values,
                )
                .unwrap();
                for (&actual, form) in values.iter().zip(&flow.integrals) {
                    assert_eq!(actual.to_bits(), quadratic(&dense(form, flow.d), &z).to_bits());
                }
            }
        }
        let mut almost_linear = dense(&flow.integrals[1], flow.d);
        almost_linear[0] = f64::from_bits(1);
        assert!(!is_affine_integral(&almost_linear, flow.d));
    }

    #[test]
    fn compact_batch_preserves_first_quadratic_transition_and_dense_affine_maps() {
        let n = 32;
        let d = n + 1;
        let rate = 0.21875;
        let mut q = vec![0.0; n * n];
        for i in 0..n {
            q[i * n + i] = 1.0;
        }
        let system = PortHamiltonian::new(
            n, 0, vec![0.0; n * n], vec![rate / n as f64; n * n], Vec::new(),
            Box::new(crate::QuadraticStorage::new(q, n).unwrap()),
        ).unwrap();
        let split = super::super::split_quadratic(
            &system, LinearSplitBudget { max_modes: n, max_omega_dt: 2.8 },
        ).unwrap();
        let flow = PreparedQuadraticFlow::new(
            &system, &split, 1.0, &[vec![0.0; n * n]], &[vec![1.0; n]],
            DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX,
                max_multiply_terms: u64::MAX,
            },
        ).unwrap();
        let mut generator = vec![0.0; d * d];
        let mut linear_form = vec![0.0; d * d];
        for i in 0..n {
            for j in 0..n {
                generator[i * d + j] = -rate / n as f64;
            }
            linear_form[i * d + n] = 0.5;
            linear_form[n * d + i] = 0.5;
        }
        let quadratic_reference = gramian::flow_integral(
            &generator, &vec![0.0; d * d], d, 1.0,
        ).unwrap();
        let linear_reference = gramian::flow_integral(
            &generator, &linear_form, d, 1.0,
        ).unwrap();
        let bits = |values: &[f64]| values.iter().map(|value| value.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&flow.transition), bits(&quadratic_reference.transition));
        assert_ne!(bits(&flow.transition), bits(&linear_reference.transition));
        assert!(matches!(flow.integrals[0], CachedIntegral::Dense { affine: true, .. }));
        match &flow.integrals[1] {
            CachedIntegral::Linear(row) => {
                assert_eq!(row.len(), d);
                assert_eq!(bits(row), bits(&linear_reference.loss_integral[n * d..]));
            }
            CachedIntegral::Dense { .. } => panic!("linear observation must use compact storage"),
        }
        let (zero_state, zero_input) = flow.linear_integral_map(0).unwrap();
        assert_eq!(zero_state, vec![0.0; n]);
        assert!(zero_input.is_empty());
        let (state_map, input_map) = flow.linear_integral_map(1).unwrap();
        assert!(input_map.is_empty());
        let mean_decay = -(-rate).exp_m1() / rate;
        for coefficient in state_map {
            assert!((coefficient - mean_decay).abs() < 32.0 * f64::EPSILON);
        }
        let initial: Vec<_> = (0..n).map(|i| (i as f64 - 15.0) / 17.0).collect();
        let step = flow.evaluate(&initial, &[]).unwrap();
        assert_eq!(step.integrals[0], 0.0);
        let expected = mean_decay * initial.iter().sum::<f64>();
        assert!((step.integrals[1] - expected).abs() < 128.0 * f64::EPSILON);
    }

    #[test]
    fn compact_preparation_keeps_reference_budgets_and_cancellation_replay() {
        let system = crate::mass_spring_damper(1.0, 9.0, 0.7).unwrap();
        let split = super::super::split_quadratic(
            &system, LinearSplitBudget { max_modes: 2, max_omega_dt: 2.8 },
        ).unwrap();
        let forms: Vec<_> = (0..24)
            .map(|i| vec![(i as f64 - 11.0) / 13.0, (7.0 - i as f64) / 19.0])
            .collect();
        let d = 4_usize;
        let budget = DenseExponentialBudget {
            max_dimension: d,
            max_dense_bytes: d * d * 8 * (2 * forms.len() + 24),
            max_multiply_terms: (d as u64).pow(3) * 224 * forms.len() as u64,
        };
        let mut polls = 0;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system, &split, 0.03, &[], &forms, budget, || {
                polls += 1;
                false
            },
        ).unwrap();
        assert_eq!(flow.integrals.len(), forms.len());
        for form in &flow.integrals {
            match form {
                CachedIntegral::Linear(row) => assert_eq!(row.len(), d),
                CachedIntegral::Dense { .. } => panic!("linear cache must stay compact"),
            }
        }
        let stop = polls / 2;
        let mut seen = 0;
        let cancelled = PreparedQuadraticFlow::new_controlled(
            &system, &split, 0.03, &[], &forms, budget, || {
                seen += 1;
                seen >= stop
            },
        );
        assert!(matches!(cancelled, Err(PreparedStepError::Cancelled)));
        let replay = PreparedQuadraticFlow::new(
            &system, &split, 0.03, &[], &forms, budget,
        ).unwrap();
        let first = flow.evaluate(&[0.3, -0.2], &[0.4]).unwrap();
        let second = replay.evaluate(&[0.3, -0.2], &[0.4]).unwrap();
        assert_eq!(first.state, second.state);
        assert_eq!(first.integrals, second.integrals);
        let mut too_small = budget;
        too_small.max_dense_bytes -= 1;
        assert!(PreparedQuadraticFlow::new(
            &system, &split, 0.03, &[], &forms, too_small,
        ).is_err());
        too_small = budget;
        too_small.max_multiply_terms -= 1;
        assert!(PreparedQuadraticFlow::new(
            &system, &split, 0.03, &[], &forms, too_small,
        ).is_err());
    }

    #[test]
    fn tiny_dimension_many_outputs_fit_the_unchanged_reference_byte_budget() {
        let system = PortHamiltonian::new(
            1, 0, vec![0.0], vec![1.0], Vec::new(),
            Box::new(crate::QuadraticStorage::new(vec![1.0], 1).unwrap()),
        ).unwrap();
        let split = super::super::split_quadratic(
            &system, LinearSplitBudget { max_modes: 1, max_omega_dt: 2.8 },
        ).unwrap();
        let forms: Vec<_> = (0..128).map(|i| vec![(i as f64 - 63.0) / 17.0]).collect();
        let d = 2_usize;
        let flow = PreparedQuadraticFlow::new(
            &system, &split, 0.1, &[], &forms,
            DenseExponentialBudget {
                max_dimension: d,
                max_dense_bytes: d * d * 8 * (2 * forms.len() + 24),
                max_multiply_terms: (d as u64).pow(3) * 224 * forms.len() as u64,
            },
        ).unwrap();
        let step = flow.evaluate(&[0.25], &[]).unwrap();
        assert!((step.state[0] - 0.25 * (-0.1_f64).exp()).abs() < 8.0 * f64::EPSILON);
        for (index, form) in flow.integrals.iter().enumerate() {
            match form {
                CachedIntegral::Linear(row) => assert_eq!(row.len(), d),
                CachedIntegral::Dense { .. } => panic!("tiny fallback must still compact its cache"),
            }
            let expected = forms[index][0] * 0.25 * -(-0.1_f64).exp_m1();
            assert!((step.integrals[index] - expected).abs() < 32.0 * f64::EPSILON);
        }
    }

}
