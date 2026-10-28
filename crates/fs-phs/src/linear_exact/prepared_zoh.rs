//! Complete held-input dynamics cached through the existing quadratic owner.
use super::{DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    finite, prepared::PreparedQuadraticFlow, product, refuse, split_quadratic};
use crate::{PhsError, PortHamiltonian, StepRecord, storage_energy_change};

/// Conservative numerical energy scale, invariant under adding a constant to H.
/// This coordinate-dependent norm is not itself thermodynamic stored energy.
pub(crate) fn energetic_scale(system: &PortHamiltonian, state: &[f64]) -> f64 {
    let mut effort = vec![0.0; state.len()];
    system.storage.gradient(state, &mut effort);
    let mut scale = 0.0;
    for (&x, &e) in state.iter().zip(&effort) {
        let term = (x * e).abs();
        if !term.is_finite() { return f64::INFINITY; }
        scale += term;
    }
    scale
}

/// Prepared complete quadratic flow, with independent heat and port integrals.
/// Preparation performs dense exponentials once. Warm steps use matrix-vector
/// products and quadratic forms; returned vectors allocate. No modes are removed.
pub struct PreparedExactZoh<'a> {
    system: &'a PortHamiltonian,
    duration: f64,
    flow: Option<PreparedQuadraticFlow>,
}

impl<'a> PreparedExactZoh<'a> {
    /// Prepare the original quadratic system for one fixed nonnegative duration.
    /// The original decomposition and complete augmented dimension caps apply.
    pub fn new(system: &'a PortHamiltonian, duration: f64,
        budget: LinearSplitBudget) -> Result<Self, PhsError> {
        if !duration.is_finite() || duration < 0.0 {
            return Err(refuse("prepared ZOH requires a finite nonnegative duration"));
        }
        let split = split_quadratic(system, budget)?;
        if !split.quadratic.exact {
            return Err(refuse("prepared ZOH requires identically quadratic storage"));
        }
        if duration == 0.0 { return Ok(Self { system, duration, flow: None }); }
        let (n, m) = (system.n, system.m);
        let q = split.quadratic_matrix();
        let loss = product(q, &product(&system.r, q, n), n);
        let outputs: Vec<Vec<f64>> = (0..m).map(|p| (0..n).map(|j|
            (0..n).map(|k| system.g[k*m+p]*q[k*n+j]).sum()).collect()).collect();
        let flow = PreparedQuadraticFlow::new(system, &split, duration, &[loss], &outputs,
            DenseExponentialBudget { max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX })?;
        Ok(Self { system, duration, flow: Some(flow) })
    }

    /// Evolve one held input. Heat is the original integral of grad(H)^T R grad(H),
    /// and output is its actual trajectory average, not the final port effort.
    pub fn step(&self, initial: &[f64], input: &[f64]) -> Result<StepRecord, PhsError> {
        let sys = self.system;
        if initial.len() != sys.n || input.len() != sys.m || !finite(initial) || !finite(input) {
            return Err(refuse("prepared ZOH requires complete finite state and input"));
        }
        let Some(flow) = &self.flow else {
            let y = sys.output(initial);
            if !finite(&y) || !sys.hamiltonian(initial).is_finite() {
                return Err(refuse("prepared ZOH initial energy/output is nonfinite"));
            }
            return Ok(StepRecord { x: initial.to_vec(), y, delta_h: 0.0,
                dissipated: 0.0, supplied: 0.0, newton_iters: 0, solver_residual: 0.0,
                solver_certificate: None });
        };
        let result = flow.evaluate(initial, input)?;
        let y: Vec<f64> = result.integrals[1..].iter().map(|v| v/self.duration).collect();
        let supplied: f64 = input.iter().zip(&result.integrals[1..]).map(|(u, integral)| u*integral).sum();
        let dissipated = result.integrals[0];
        let delta_h = storage_energy_change(sys.storage.as_ref(), initial, &result.state,
            &mut vec![0.0; sys.n], &mut vec![0.0; sys.n]);
        let scale = energetic_scale(sys, initial).max(energetic_scale(sys, &result.state))
            .max(delta_h.abs()).max(supplied.abs()).max(dissipated.abs());
        let residual = delta_h + dissipated - supplied;
        if !finite(&y) || !sys.hamiltonian(&result.state).is_finite()
            || !finite(&[delta_h, supplied, dissipated, scale, residual])
            || dissipated < -1e-12*scale || residual.abs() > 1e-12*scale {
            return Err(refuse("prepared ZOH physical energy balance refused"));
        }
        Ok(StepRecord { x: result.state, y, delta_h, dissipated, supplied,
            newton_iters: 0, solver_residual: 0.0, solver_certificate: None })
    }
}
