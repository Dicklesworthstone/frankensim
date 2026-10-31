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
    CommonFlowProjection, DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow,
    prepared_zoh::energetic_scale_into, refuse, split_quadratic,
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
    projection: CommonFlowProjection,
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
            resistance_witness,
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
            resistance_witness,
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
        let projection = CommonFlowProjection::new(
            &original, &b, 1, budget.max_modes, &mut || false,
        ).map_err(|error| match error {
            crate::PreparedStepError::Solver(error) => error,
            crate::PreparedStepError::Cancelled => refuse("uncancelled common-flow preparation interrupted"),
        })?;
        let mut tangent = original.share_storage();
        let (j, r, g) = projection.structure();
        tangent.j.copy_from_slice(j);
        tangent.r.copy_from_slice(r);
        tangent.g.copy_from_slice(g);
        // The original algebraic multiplier is reconstructed at the endpoint;
        // its zero-storage tangent coordinate has exactly zero evolution.
        for i in 0..n {
            tangent.j[i * n + a] = 0.0;
            tangent.j[a * n + i] = 0.0;
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
            resistance_witness,
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
                    resistance_witness,
                },
                tangent,
                flow,
                heat,
                projection,
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
        self.core
            .admit(state, input, &mut vec![0.0; self.state_dim()])?;
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
    fn constraint(&self, state: &[f64], _effort: &mut [f64]) -> Result<(), PhsError> {
        self.projection.admit_state(state)
    }
    fn admit(&self, state: &[f64], input: &[f64], effort: &mut [f64]) -> Result<(), PhsError> {
        if state.len() != self.original.n
            || input.len() != self.original.m
            || !finite(state)
            || !finite(input)
        {
            return Err(refuse(
                "common-flow requires complete finite state and input",
            ));
        }
        self.constraint(state, effort)
    }
    fn reaction(&self, state: &[f64], input: &[f64]) -> Result<f64, PhsError> {
        let mut value = [0.0];
        self.projection.reaction_into(state, input, &mut value)?;
        Ok(value[0])
    }
    fn evaluate(
        &self,
        initial: &[f64],
        input: &[f64],
        scratch: &mut Scratch,
    ) -> Result<PreparedStepRecord, PhsError> {
        self.admit(initial, input, &mut scratch.effort)?;
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
        self.constraint(&scratch.endpoint, &mut scratch.effort)?;
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
