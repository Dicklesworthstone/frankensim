//! Complete quadratic evolution with one pure common-flow constraint.
//!
//! For differential state x and multiplier lambda, the original equations are
//! x'=(J-R)Qx+B lambda+Gu, B^T Qx=0. With S=B^T QB>0 and
//! P=I-B S^-1 B^T Q, the tangent dynamics are P(J-R)P^T Qx+PGu.
//! On the admitted constraint P^T Qx=Qx, so this is the original DAE, without
//! dropping a coordinate or projecting a physical state. Original G observes
//! port work and original R produces heat; lambda is the differentiated
//! constraint reaction, not a stored-energy coordinate.
use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow,
    prepared_zoh::energetic_scale_into, product, refuse, split_quadratic,
};
use crate::{
    DescriptorPortHamiltonian, PhsError, PortHamiltonian, PreparedStepRecord, StepRecord,
    storage_energy_change,
};

struct Scratch {
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
    output: Vec<f64>,
    midpoint: Vec<f64>,
    effort: Vec<f64>,
    heat_coordinates: Vec<f64>,
}
impl Scratch {
    fn new(n: usize, m: usize) -> Self {
        Self {
            augmented: vec![0.0; n + m + 1],
            balanced: vec![0.0; n],
            endpoint: vec![0.0; n],
            integrals: vec![0.0; m],
            output: vec![0.0; m],
            midpoint: vec![0.0; n],
            effort: vec![0.0; n],
            heat_coordinates: vec![0.0; 2 * n],
        }
    }
}

struct Core {
    original: DescriptorPortHamiltonian,
    tangent: PortHamiltonian,
    flow: Option<PreparedQuadraticFlow>,
    heat: Option<PositiveHeat>,
    constraint: Vec<f64>,
    reaction_state: Vec<f64>,
    reaction_input: Vec<f64>,
    dt: f64,
}

/// Prepared exact held-input flow of one quadratic common-velocity junction.
///
/// All original physical coordinates, ports and material histories remain in
/// their original order, followed by the original single multiplier. A finite
/// incoming multiplier carries no energy and is replaced by the instantaneous
/// ENDPOINT reaction for the held input. Initial differential state must already
/// satisfy the original constraint. No projection or drift repair occurs.
///
/// Preparation reuses the shared complete quadratic propagator and positive
/// physical-heat integral. `step_into` reuses scratch without owner allocation
/// or exponentials. This is not a general descriptor solver: nonlinear storage,
/// multiple constraints, algebraic resistance/sources and a singular effective
/// inverse inertia are refused. Consumer source bandwidth remains explicit.
pub struct PreparedCommonFlowZoh {
    core: Core,
    scratch: Scratch,
}

impl PreparedCommonFlowZoh {
    /// Own the original single-constraint descriptor and one fixed-clock cache.
    /// The complete state, held-input and dense preparation limits apply before
    /// copying matrices; no modal or material-memory truncation is permitted.
    pub fn new_owned(
        system: DescriptorPortHamiltonian,
        duration: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        let DescriptorPortHamiltonian {
            n,
            n_diff,
            m,
            j,
            r,
            g,
            storage,
        } = system;
        if n_diff == 0
            || n_diff.checked_add(1) != Some(n)
            || !duration.is_finite()
            || duration < 0.0
        {
            return Err(refuse(
                "common-flow ZOH requires one constraint and a finite nonnegative clock",
            ));
        }
        let mut original = PortHamiltonian {
            n,
            m,
            j,
            r,
            g,
            storage,
        };
        let original_split = split_quadratic(&original, budget)?;
        let q = original_split.quadratic_matrix();
        let a = n_diff;
        if !original_split.is_exact_quadratic()
            || (0..n).any(|i| {
                q[a * n + i] != 0.0
                    || q[i * n + a] != 0.0
                    || original.r[a * n + i] != 0.0
                    || original.r[i * n + a] != 0.0
            })
            || original.j[a * n + a] != 0.0
            || (0..m).any(|p| original.g[a * m + p] != 0.0)
        {
            return Err(refuse(
                "common-flow ZOH requires exact quadratic differential storage and a pure algebraic constraint",
            ));
        }
        let b: Vec<f64> = (0..n).map(|i| original.j[i * n + a]).collect();
        let constraint: Vec<f64> = (0..n)
            .map(|j| (0..n).map(|i| b[i] * q[i * n + j]).sum())
            .collect();
        let inertia: f64 = b.iter().zip(&constraint).map(|(b, c)| b * c).sum();
        if !finite(&constraint) || !inertia.is_finite() || inertia <= 0.0 {
            return Err(refuse(
                "common-flow constraint requires positive representable effective inverse inertia",
            ));
        }
        let normalized: Vec<f64> = constraint.iter().map(|c| c / inertia).collect();
        if !finite(&normalized) {
            return Err(refuse("common-flow constraint scaling is unrepresentable"));
        }
        // The original scalar differentiated constraint determines the reaction.
        // These are exact declared linear maps, never response probes.
        let reaction_effort: Vec<f64> = (0..n)
            .map(|j| {
                -(0..n_diff)
                    .map(|i| normalized[i] * (original.j[i * n + j] - original.r[i * n + j]))
                    .sum::<f64>()
            })
            .collect();
        let reaction_state: Vec<f64> = (0..n)
            .map(|j| (0..n).map(|i| reaction_effort[i] * q[i * n + j]).sum())
            .collect();
        let reaction_input: Vec<f64> = (0..m)
            .map(|p| {
                -(0..n_diff)
                    .map(|i| normalized[i] * original.g[i * m + p])
                    .sum::<f64>()
            })
            .collect();
        let mut projector = vec![0.0; n * n];
        for i in 0..n_diff {
            for j in 0..n_diff {
                projector[i * n + j] = if i == j { 1.0 } else { 0.0 } - b[i] * normalized[j];
            }
        }
        // The multiplier has no evolution or storage in the tangent flow. Its
        // returned value is reconstructed after the complete physical step.
        projector[a * n + a] = 1.0;
        let mut transpose = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                transpose[i * n + j] = projector[j * n + i];
            }
        }
        let mut tangent = original.share_storage();
        for i in 0..n {
            tangent.j[i * n + a] = 0.0;
            tangent.j[a * n + i] = 0.0;
        }
        tangent.j = product(&product(&projector, &tangent.j, n), &transpose, n);
        tangent.r = product(&product(&projector, &original.r, n), &transpose, n);
        // Stamp each analytical skew/symmetric pair from the same contraction,
        // avoiding unequal roundoff in the two evaluations of one coefficient.
        for i in 0..n {
            tangent.j[i * n + i] = 0.0;
            for j in i + 1..n {
                tangent.j[j * n + i] = -tangent.j[i * n + j];
                tangent.r[j * n + i] = tangent.r[i * n + j];
            }
            for p in 0..m {
                tangent.g[i * m + p] = (0..n)
                    .map(|k| projector[i * n + k] * original.g[k * m + p])
                    .sum();
            }
        }
        if !finite(&reaction_state) || !finite(&reaction_input) || !finite(&projector) {
            return Err(refuse("common-flow tangent or reaction is unrepresentable"));
        }
        let split = split_quadratic(&tangent, budget)?;
        let (flow, heat) = if duration == 0.0 {
            (None, None)
        } else {
            let outputs: Vec<Vec<f64>> = (0..m)
                .map(|p| {
                    (0..n)
                        .map(|j| (0..n).map(|i| original.g[i * m + p] * q[i * n + j]).sum())
                        .collect()
                })
                .collect();
            let flow = PreparedQuadraticFlow::new(
                &tangent,
                &split,
                duration,
                &[],
                &outputs,
                DenseExponentialBudget {
                    max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                    max_dense_bytes: usize::MAX,
                    max_multiply_terms: u64::MAX,
                },
            )?;
            // Heat belongs to ORIGINAL material R, along the constrained flow.
            // Factoring the projected singular resistance instead would both
            // lose that physical meaning and introduce artificial null modes.
            let heat = PositiveHeat::new(&tangent, &split, duration, &original.r)?;
            (Some(flow), Some(heat))
        };
        let PortHamiltonian {
            n,
            m,
            j,
            r,
            g,
            storage,
        } = original;
        Ok(Self {
            core: Core {
                original: DescriptorPortHamiltonian {
                    n,
                    n_diff,
                    m,
                    j,
                    r,
                    g,
                    storage,
                },
                tangent,
                flow,
                heat,
                constraint,
                reaction_state,
                reaction_input,
                dt: duration,
            },
            scratch: Scratch::new(n, m),
        })
    }

    /// The unchanged original physical descriptor and its constitutive law.
    pub fn system(&self) -> &DescriptorPortHamiltonian {
        &self.core.original
    }
    pub fn state_dim(&self) -> usize {
        self.core.original.n
    }
    pub fn differential_dim(&self) -> usize {
        self.core.original.n_diff
    }
    pub fn port_dim(&self) -> usize {
        self.core.original.m
    }
    pub fn duration(&self) -> f64 {
        self.core.dt
    }
    pub fn hamiltonian(&self, state: &[f64]) -> f64 {
        self.core.original.hamiltonian(state)
    }

    /// Original multiplier at this physical state and held external input.
    /// A finite incoming multiplier does not affect the differentiated constraint.
    pub fn instantaneous_reaction(&self, state: &[f64], input: &[f64]) -> Result<f64, PhsError> {
        self.core.admit(state, input)?;
        self.core.reaction(state, input)
    }

    /// Allocate the complete candidate and original mean port outputs.
    pub fn step(&self, initial: &[f64], input: &[f64]) -> Result<StepRecord, PhsError> {
        let mut scratch = Scratch::new(self.state_dim(), self.port_dim());
        let record = self.core.evaluate(initial, input, &mut scratch)?;
        Ok(StepRecord {
            x: scratch.endpoint,
            y: scratch.output,
            delta_h: record.delta_h,
            dissipated: record.dissipated,
            supplied: record.supplied,
            newton_iters: record.newton_iters,
            solver_residual: record.solver_residual,
            solver_certificate: record.solver_certificate,
        })
    }

    /// Same physical trajectory and independent accounts using existing scratch.
    /// Both caller output buffers remain unchanged on every refusal.
    pub fn step_into(
        &mut self,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        if next.len() != self.state_dim() || average_output.len() != self.port_dim() {
            return Err(refuse(
                "common-flow output buffers must match the complete state and ports",
            ));
        }
        let record = self.core.evaluate(initial, input, &mut self.scratch)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok(record)
    }
}

impl Core {
    fn constraint(&self, state: &[f64]) -> Result<(), PhsError> {
        let mut residual = 0.0;
        let mut scale = 0.0;
        for (c, x) in self.constraint.iter().zip(state) {
            let term = c * x;
            residual += term;
            scale += term.abs();
        }
        if !finite(&[residual, scale])
            || residual.abs() > 64.0 * f64::EPSILON * self.original.n as f64 * scale
        {
            return Err(refuse(
                "common-flow original velocity constraint is incompatible",
            ));
        }
        Ok(())
    }
    fn admit(&self, state: &[f64], input: &[f64]) -> Result<(), PhsError> {
        if state.len() != self.original.n
            || input.len() != self.original.m
            || !finite(state)
            || !finite(input)
        {
            return Err(refuse(
                "common-flow requires complete finite state and input",
            ));
        }
        self.constraint(state)
    }
    fn reaction(&self, state: &[f64], input: &[f64]) -> Result<f64, PhsError> {
        let value = self
            .reaction_state
            .iter()
            .zip(state)
            .map(|(a, x)| a * x)
            .sum::<f64>()
            + self
                .reaction_input
                .iter()
                .zip(input)
                .map(|(a, u)| a * u)
                .sum::<f64>();
        if !value.is_finite() {
            return Err(refuse("common-flow instantaneous reaction is nonfinite"));
        }
        Ok(value)
    }
    fn evaluate(
        &self,
        initial: &[f64],
        input: &[f64],
        scratch: &mut Scratch,
    ) -> Result<PreparedStepRecord, PhsError> {
        self.admit(initial, input)?;
        let (n, m) = (self.original.n, self.original.m);
        let (delta_h, dissipated, supplied) = if let Some(flow) = &self.flow {
            flow.evaluate_into(
                initial,
                input,
                &mut scratch.augmented,
                &mut scratch.balanced,
                &mut scratch.endpoint,
                &mut scratch.integrals,
            )?;
            for (out, integral) in scratch.output.iter_mut().zip(&scratch.integrals) {
                *out = integral / self.dt;
            }
            let supplied: f64 = input
                .iter()
                .zip(&scratch.integrals)
                .map(|(u, y)| u * y)
                .sum();
            let dissipated = self
                .heat
                .as_ref()
                .ok_or_else(|| refuse("common-flow physical heat cache missing"))?
                .evaluate_into(
                    &self.tangent,
                    initial,
                    input,
                    &mut scratch.effort,
                    &mut scratch.balanced,
                    &mut scratch.heat_coordinates,
                )?;
            let delta_h = storage_energy_change(
                self.original.storage.as_ref(),
                initial,
                &scratch.endpoint,
                &mut scratch.midpoint,
                &mut scratch.effort,
            );
            (delta_h, dissipated, supplied)
        } else {
            scratch.endpoint.copy_from_slice(initial);
            self.original.storage.gradient(initial, &mut scratch.effort);
            for p in 0..m {
                scratch.output[p] = (0..n)
                    .map(|i| self.original.g[i * m + p] * scratch.effort[i])
                    .sum();
            }
            (0.0, 0.0, 0.0)
        };
        self.constraint(&scratch.endpoint)?;
        scratch.endpoint[self.original.n_diff] = self.reaction(&scratch.endpoint, input)?;
        let scale = energetic_scale_into(&self.tangent, initial, &mut scratch.effort)
            .max(energetic_scale_into(
                &self.tangent,
                &scratch.endpoint,
                &mut scratch.effort,
            ))
            .max(delta_h.abs())
            .max(supplied.abs())
            .max(dissipated.abs());
        let residual = delta_h + dissipated - supplied;
        if !finite(&scratch.endpoint)
            || !finite(&scratch.output)
            || !self.original.hamiltonian(initial).is_finite()
            || !self.original.hamiltonian(&scratch.endpoint).is_finite()
            || !finite(&[delta_h, dissipated, supplied, scale, residual])
            || dissipated < -1e-12 * scale
            || residual.abs() > 1e-12 * scale
        {
            return Err(refuse(
                "common-flow original physical energy balance refused",
            ));
        }
        Ok(PreparedStepRecord {
            delta_h,
            dissipated,
            supplied,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        })
    }
}
