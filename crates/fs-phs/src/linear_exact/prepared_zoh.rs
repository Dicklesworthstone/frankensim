//! Complete held-input dynamics cached through the existing quadratic owner.
use super::{DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    finite, prepared::PreparedQuadraticFlow, product, refuse, split_quadratic};
use crate::{PhsError, PortHamiltonian, PreparedStepRecord, StepRecord, storage_energy_change};

/// Conservative numerical energy scale, invariant under adding a constant to H.
/// This coordinate-dependent norm is not itself thermodynamic stored energy.
pub(crate) fn energetic_scale(system: &PortHamiltonian, state: &[f64]) -> f64 {
    energetic_scale_into(system, state, &mut vec![0.0; state.len()])
}
fn energetic_scale_into(system: &PortHamiltonian, state: &[f64], effort: &mut [f64]) -> f64 {
    effort.fill(0.0);
    system.storage.gradient(state, effort);
    let mut scale = 0.0;
    for (&x, &e) in state.iter().zip(effort.iter()) {
        let term = (x * e).abs();
        if !term.is_finite() { return f64::INFINITY; }
        scale += term;
    }
    scale
}

enum System<'a> { Borrowed(&'a PortHamiltonian), Owned(Box<PortHamiltonian>) }
impl System<'_> {
    fn get(&self) -> &PortHamiltonian {
        match self { Self::Borrowed(system) => system, Self::Owned(system) => system }
    }
}
struct Scratch {
    augmented: Vec<f64>, balanced: Vec<f64>, endpoint: Vec<f64>, integrals: Vec<f64>,
    midpoint: Vec<f64>, effort: Vec<f64>, output: Vec<f64>,
}
impl Scratch {
    fn new(n: usize, m: usize) -> Self {
        Self { augmented: vec![0.0; n+m+1], balanced: vec![0.0; n], endpoint: vec![0.0; n],
            integrals: vec![0.0; m+1], midpoint: vec![0.0; n], effort: vec![0.0; n],
            output: vec![0.0; m] }
    }
}

/// Prepared complete quadratic flow, with independent heat and port integrals.
/// Preparation performs dense exponentials once. Warm steps use matrix-vector
/// products and quadratic forms. `step` returns allocating vectors; `step_into`
/// reuses complete scratch with no owner allocations or matrix exponentials.
pub struct PreparedExactZoh<'a> {
    system: System<'a>,
    duration: f64,
    flow: Option<PreparedQuadraticFlow>,
    scratch: Scratch,
}

impl PreparedExactZoh<'static> {
    /// Own the original immutable system and its fixed-clock cache together.
    /// This permits a movable render runtime without self-referential borrowing,
    /// copying storage, changing its state order, or rebuilding per sample.
    pub fn new_owned(system: PortHamiltonian, duration: f64,
        budget: LinearSplitBudget) -> Result<Self, PhsError> {
        Self::prepare(System::Owned(Box::new(system)), duration, budget)
    }
}
impl<'a> PreparedExactZoh<'a> {
    /// Prepare the original quadratic system for one fixed nonnegative duration.
    /// The original decomposition and complete augmented dimension caps apply.
    pub fn new(system: &'a PortHamiltonian, duration: f64,
        budget: LinearSplitBudget) -> Result<Self, PhsError> {
        Self::prepare(System::Borrowed(system), duration, budget)
    }
    fn prepare(owner: System<'a>, duration: f64, budget: LinearSplitBudget) -> Result<Self, PhsError> {
        let system = owner.get();
        if !duration.is_finite() || duration < 0.0 {
            return Err(refuse("prepared ZOH requires a finite nonnegative duration"));
        }
        let split = split_quadratic(system, budget)?;
        if !split.quadratic.exact {
            return Err(refuse("prepared ZOH requires identically quadratic storage"));
        }
        let (n, m) = (system.n, system.m);
        let scratch = Scratch::new(n, m);
        if duration == 0.0 { return Ok(Self { system: owner, duration, flow: None, scratch }); }
        let q = split.quadratic_matrix();
        let loss = product(q, &product(&system.r, q, n), n);
        let outputs: Vec<Vec<f64>> = (0..m).map(|p| (0..n).map(|j|
            (0..n).map(|k| system.g[k*m+p]*q[k*n+j]).sum()).collect()).collect();
        let flow = PreparedQuadraticFlow::new(system, &split, duration, &[loss], &outputs,
            DenseExponentialBudget { max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX })?;
        Ok(Self { system: owner, duration, flow: Some(flow), scratch })
    }

    /// Original immutable physical owner, including complete material storage.
    pub fn system(&self) -> &PortHamiltonian { self.system.get() }
    /// Prepared held-input duration [s]. It never changes during a warm step.
    pub fn duration(&self) -> f64 { self.duration }

    /// Evolve one held input. Heat is the original integral of grad(H)^T R grad(H),
    /// and output is its actual trajectory average, not the final port effort.
    pub fn step(&self, initial: &[f64], input: &[f64]) -> Result<StepRecord, PhsError> {
        let sys = self.system();
        let mut scratch = Scratch::new(sys.n, sys.m);
        let record = evaluate_step(sys, self.duration, self.flow.as_ref(), initial, input, &mut scratch)?;
        Ok(StepRecord { x: scratch.endpoint, y: scratch.output, delta_h: record.delta_h,
            dissipated: record.dissipated, supplied: record.supplied,
            newton_iters: record.newton_iters, solver_residual: record.solver_residual,
            solver_certificate: record.solver_certificate })
    }

    /// Same exact transition and independently integrated heat/work, into caller
    /// buffers. No owner allocation or exponential occurs on a successful warm
    /// call; storage callbacks must obey the existing allocation-free contract.
    /// Neither output is changed on ANY refusal. Scratch may change.
    pub fn step_into(&mut self, initial: &[f64], input: &[f64], next: &mut [f64],
        average_output: &mut [f64]) -> Result<PreparedStepRecord, PhsError> {
        let sys = self.system.get();
        if next.len() != sys.n || average_output.len() != sys.m {
            return Err(refuse("prepared ZOH output buffers must match the complete state and ports"));
        }
        let record = evaluate_step(sys, self.duration, self.flow.as_ref(), initial, input, &mut self.scratch)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok(record)
    }
}

fn evaluate_step(sys: &PortHamiltonian, duration: f64, flow: Option<&PreparedQuadraticFlow>,
    initial: &[f64], input: &[f64], scratch: &mut Scratch) -> Result<PreparedStepRecord, PhsError> {
    if initial.len() != sys.n || input.len() != sys.m || !finite(initial) || !finite(input) {
        return Err(refuse("prepared ZOH requires complete finite state and input"));
    }
    let Some(flow) = flow else {
        scratch.effort.fill(0.0);
        sys.storage.gradient(initial, &mut scratch.effort);
        for p in 0..sys.m {
            scratch.output[p] = (0..sys.n).map(|i| sys.g[i*sys.m+p]*scratch.effort[i]).sum();
        }
        if !finite(&scratch.output) || !sys.hamiltonian(initial).is_finite() {
            return Err(refuse("prepared ZOH initial energy/output is nonfinite"));
        }
        scratch.endpoint.copy_from_slice(initial);
        return Ok(PreparedStepRecord { delta_h: 0.0, dissipated: 0.0, supplied: 0.0,
            newton_iters: 0, solver_residual: 0.0, solver_certificate: None });
    };
    flow.evaluate_into(initial, input, &mut scratch.augmented, &mut scratch.balanced,
        &mut scratch.endpoint, &mut scratch.integrals)?;
    for (out, value) in scratch.output.iter_mut().zip(&scratch.integrals[1..]) { *out = value/duration; }
    let supplied: f64 = input.iter().zip(&scratch.integrals[1..]).map(|(u, integral)| u*integral).sum();
    let dissipated = scratch.integrals[0];
    let delta_h = storage_energy_change(sys.storage.as_ref(), initial, &scratch.endpoint,
        &mut scratch.midpoint, &mut scratch.effort);
    let scale = energetic_scale_into(sys, initial, &mut scratch.effort)
        .max(energetic_scale_into(sys, &scratch.endpoint, &mut scratch.effort))
        .max(delta_h.abs()).max(supplied.abs()).max(dissipated.abs());
    let residual = delta_h + dissipated - supplied;
    if !finite(&scratch.output) || !sys.hamiltonian(&scratch.endpoint).is_finite()
        || !finite(&[delta_h, supplied, dissipated, scale, residual])
        || dissipated < -1e-12*scale || residual.abs() > 1e-12*scale {
        return Err(refuse("prepared ZOH physical energy balance refused"));
    }
    Ok(PreparedStepRecord { delta_h, dissipated, supplied,
        newton_iters: 0, solver_residual: 0.0, solver_certificate: None })
}
