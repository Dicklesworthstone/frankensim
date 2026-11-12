//! Complete held-input dynamics cached through the existing quadratic owner.
use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, balance_admitted,
    finite, positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow, refuse, split_quadratic,
};
use crate::interconnection::accumulation::Sum;
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, StepRecord,
    storage_energy_change,
};

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

fn uncancelled(error: PreparedStepError) -> PhsError {
    match error {
        PreparedStepError::Solver(error) => error,
        PreparedStepError::Cancelled => refuse("uncancelled exact ZOH operation interrupted"),
    }
}

/// Conservative numerical energy scale, invariant under adding a constant to H.
/// This coordinate-dependent norm is not itself thermodynamic stored energy.
pub(crate) fn energetic_scale(system: &PortHamiltonian, state: &[f64]) -> f64 {
    energetic_scale_into(system, state, &mut vec![0.0; state.len()])
}
pub(super) fn energetic_scale_into(system: &PortHamiltonian, state: &[f64], effort: &mut [f64]) -> f64 {
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

pub(super) enum System<'a> { Borrowed(&'a PortHamiltonian), Owned(Box<PortHamiltonian>) }
impl System<'_> {
    pub(super) fn get(&self) -> &PortHamiltonian {
        match self { Self::Borrowed(system) => system, Self::Owned(system) => system }
    }
}
#[derive(Clone)]
struct Scratch {
    augmented: Vec<f64>, balanced: Vec<f64>, endpoint: Vec<f64>, integrals: Vec<f64>,
    midpoint: Vec<f64>, effort: Vec<f64>, output: Vec<f64>,
    heat_state: Vec<f64>,
    port_integrals: Vec<Sum>,
}
impl Scratch {
    fn new(n: usize, m: usize) -> Self {
        Self { augmented: vec![0.0; n+m+1], balanced: vec![0.0; n], endpoint: vec![0.0; n],
            integrals: vec![0.0; m], midpoint: vec![0.0; n], effort: vec![0.0; n],
            output: vec![0.0; m], heat_state: vec![0.0; 2*n],
            port_integrals: vec![Sum::default(); m] }
    }
}

/// Reusable caller-owned scratch for one complete state/port layout.
///
/// This contains no accepted physical state. Clones have independent buffers,
/// and one workspace can be reused with any immutable [`PreparedExactZoh`]
/// having the same state and port dimensions, even when its operator or clock
/// differs. A refused evaluation may change scratch but never caller outputs.
#[derive(Clone)]
pub struct ExactZohWorkspace {
    scratch: Scratch,
}

impl core::fmt::Debug for ExactZohWorkspace {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ExactZohWorkspace")
            .field("state_dim", &self.scratch.endpoint.len())
            .field("port_dim", &self.scratch.output.len())
            .finish_non_exhaustive()
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
    heat: Option<PositiveHeat>,
    // Original-coordinate maps for each interval-mean output. These are
    // extracted from the same cached integral, never from trial trajectories.
    mean_state: Vec<f64>,
    mean_input: Vec<f64>,
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

    /// Own the same complete fixed-duration cache with cooperative cancellation.
    /// The original split/augmented dimension and work admissions remain. Polls
    /// bracket native splitting and positive-heat preparation and run inside the
    /// quadratic-flow preparation. Splitting, storage callbacks and one original
    /// positive-heat construction remain bounded, indivisible operations.
    pub fn new_owned_controlled(
        system: PortHamiltonian,
        duration: f64,
        budget: LinearSplitBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare_controlled(System::Owned(Box::new(system)), duration, budget, cancelled)
    }
}
impl<'a> PreparedExactZoh<'a> {
    /// Prepare the original quadratic system for one fixed nonnegative duration.
    /// The original decomposition and complete augmented dimension caps apply.
    pub fn new(system: &'a PortHamiltonian, duration: f64,
        budget: LinearSplitBudget) -> Result<Self, PhsError> {
        Self::prepare(System::Borrowed(system), duration, budget)
    }

    /// Borrow the original owner with the same cancellation scope and admission
    /// as [`Self::new_owned_controlled`]. No partial cache is returned.
    pub fn new_controlled(
        system: &'a PortHamiltonian,
        duration: f64,
        budget: LinearSplitBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare_controlled(System::Borrowed(system), duration, budget, cancelled)
    }

    fn prepare(
        owner: System<'a>,
        duration: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        Self::prepare_controlled(owner, duration, budget, || false).map_err(uncancelled)
    }

    fn prepare_controlled(
        owner: System<'a>,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let system = owner.get();
        if !duration.is_finite() || duration < 0.0 {
            return Err(refuse("prepared ZOH requires a finite nonnegative duration").into());
        }
        let split = split_quadratic(system, budget)?;
        poll(&mut cancelled)?;
        if !split.quadratic.exact {
            return Err(refuse("prepared ZOH requires identically quadratic storage").into());
        }
        let (n, m) = (system.n, system.m);
        let scratch = Scratch::new(n, m);
        let q = split.quadratic_matrix();
        let mut outputs = Vec::with_capacity(m);
        for p in 0..m {
            poll(&mut cancelled)?;
            outputs.push(
                (0..n)
                    .map(|j| {
                        crate::port_output_component(&system.g, m, p, |k| q[k * n + j]).ok_or_else(
                            || refuse("prepared port observation map is unrepresentable"),
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
        if duration == 0.0 {
            let mean_state: Vec<f64> = outputs.into_iter().flatten().collect();
            if !finite(&mean_state) {
                return Err(refuse("prepared ZOH instantaneous output map is nonfinite").into());
            }
            poll(&mut cancelled)?;
            return Ok(Self { system: owner, duration, flow: None, heat: None,
                mean_state, mean_input: Vec::new(), scratch });
        }
        let flow = PreparedQuadraticFlow::new_controlled(system, &split, duration, &[], &outputs,
            DenseExponentialBudget { max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX }, &mut cancelled)?;
        let mut mean_state = Vec::with_capacity(m*n);
        let mut mean_input = Vec::with_capacity(m*m);
        for p in 0..m {
            poll(&mut cancelled)?;
            let (state, input) = flow.linear_integral_map(p)?;
            mean_state.extend(state.into_iter().map(|value| value/duration));
            mean_input.extend(input.into_iter().map(|value| value/duration));
        }
        if !finite(&mean_state) || !finite(&mean_input) {
            return Err(refuse("prepared ZOH mean output map is nonfinite").into());
        }
        poll(&mut cancelled)?;
        let heat = PositiveHeat::new(system, &split, duration, &system.r)?;
        poll(&mut cancelled)?;
        Ok(Self { system: owner, duration, flow: Some(flow), heat: Some(heat),
            mean_state, mean_input, scratch })
    }

    /// Original immutable physical owner, including complete material storage.
    pub fn system(&self) -> &PortHamiltonian { self.system.get() }
    /// Prepared held-input duration [s]. It never changes during a warm step.
    pub fn duration(&self) -> f64 { self.duration }

    /// Allocate independent scratch once for this complete state/port layout.
    /// Immutable prepared operators can then be shared across continuing
    /// trajectories without locking or copying the exponential/heat caches.
    pub fn workspace(&self) -> ExactZohWorkspace {
        let system = self.system();
        ExactZohWorkspace { scratch: Scratch::new(system.n, system.m) }
    }

    /// Apply the unchanged exact flow and physical balance gate using external
    /// scratch. This allocates nothing and evaluates no exponentials on a warm
    /// call, subject to the same storage-callback contract as [`Self::step_into`].
    ///
    /// Both state AND port dimensions must match the workspace. All caller
    /// outputs remain unchanged on every refusal; scratch is reusable afterward.
    pub fn step_with_workspace_into(
        &self,
        workspace: &mut ExactZohWorkspace,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_with_workspace_into_controlled(
            workspace,
            initial,
            input,
            next,
            average_output,
            || false,
        )
        .map_err(uncancelled)
    }

    /// The same exact flow with cancellation before and after the original
    /// matrix-vector/port-work block, positive-heat block and energy admission.
    /// Those bounded blocks and storage callbacks are not preemptible. Caller
    /// outputs publish together only after the final poll; refused scratch is
    /// reusable with any admitted clock having the same complete layout.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_workspace_into_controlled(
        &self,
        workspace: &mut ExactZohWorkspace,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        poll(&mut cancelled)?;
        let system = self.system();
        if workspace.scratch.endpoint.len() != system.n
            || workspace.scratch.output.len() != system.m
            || next.len() != system.n || average_output.len() != system.m
        {
            return Err(refuse("prepared ZOH workspace and outputs must match complete state and ports").into());
        }
        let (record, _) = evaluate_step_with_work_inputs_controlled(
            system,
            self.duration,
            self.flow.as_ref(),
            self.heat.as_ref(),
            initial,
            input,
            [],
            &mut workspace.scratch,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&workspace.scratch.endpoint);
        average_output.copy_from_slice(&workspace.scratch.output);
        Ok(record)
    }

    /// Original physical resistance heat and its derivative with respect to
    /// every held input at FIXED initial state and duration. `gradient` has one
    /// entry per original port. This differentiates the retained positive factor
    /// directly; it forms no expanded Gramian, heat difference or new trajectory.
    /// `average_output_map` supplies the separate conjugate-output derivatives.
    ///
    /// Warm work uses existing scratch and allocates nothing. Polls bracket the
    /// original heat evaluation, then occur per derivative column/factor row.
    /// Both dimensions and source support are checked; output publishes only on
    /// success. This observation is not a replacement for the step's independent
    /// physical balance admission. Zero duration/resistance gives exact zeros.
    pub fn heat_input_gradient_into_controlled(
        &self,
        workspace: &mut ExactZohWorkspace,
        initial: &[f64],
        input: &[f64],
        gradient: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        poll(&mut cancelled)?;
        let system = self.system();
        if workspace.scratch.endpoint.len() != system.n
            || workspace.scratch.output.len() != system.m
            || initial.len() != system.n
            || input.len() != system.m
            || gradient.len() != system.m
            || !finite(initial)
            || !finite(input)
            || !system.hamiltonian(initial).is_finite()
        {
            return Err(refuse(
                "prepared heat derivative requires complete finite state, ports and workspace",
            )
            .into());
        }
        let scratch = &mut workspace.scratch;
        let heat = if let Some(heat) = &self.heat {
            heat.input_gradient_into_controlled(
                system,
                initial,
                input,
                &mut scratch.effort,
                &mut scratch.balanced,
                &mut scratch.heat_state,
                &mut scratch.output,
                &mut cancelled,
            )?
        } else {
            scratch.output.fill(0.0);
            0.0
        };
        poll(&mut cancelled)?;
        gradient.copy_from_slice(&scratch.output);
        Ok(heat)
    }

    /// Magnitude of the retained original port-work factors, before opposing
    /// ports cancel: `sum |u_p| * (|integral_high_p| + |integral_low_p|)` [J].
    /// A consumer can scale arithmetic allowance for converting those native
    /// integrals to rounded mean outputs and back to wire work. This is not an
    /// integration-error bound or a replacement physical balance tolerance.
    ///
    /// Call immediately after a SUCCESSFUL step of this owner into `workspace`,
    /// with that step's same `input`, before any other owner/clock uses scratch.
    /// Layout is checked; workspace intentionally retains no owner/input history.
    pub fn port_work_magnitude(
        &self,
        workspace: &ExactZohWorkspace,
        input: &[f64],
    ) -> Result<f64, PhsError> {
        let system = self.system();
        if workspace.scratch.endpoint.len() != system.n
            || workspace.scratch.output.len() != system.m
            || input.len() != system.m
            || !finite(input)
        {
            return Err(refuse(
                "port work magnitude requires complete matching workspace and input",
            ));
        }
        let mut magnitude = 0.0;
        for (integral, &u) in workspace.scratch.port_integrals.iter().zip(input) {
            let (high, low) = integral
                .components()
                .ok_or_else(|| refuse("retained port integral is unrepresentable"))?;
            magnitude += u.abs() * high.abs();
            magnitude += u.abs() * low.abs();
        }
        if !magnitude.is_finite() {
            return Err(refuse("port work magnitude is unrepresentable"));
        }
        Ok(magnitude)
    }

    /// Original cached coefficients of one actual interval-mean output:
    /// `y_bar = state_row * initial + input_row * held_input`.
    /// Rows retain the original coordinate/port ordering. At zero duration
    /// the input row is empty and the state row is the instantaneous output.
    /// This immutable view performs no propagation, allocation or fitting;
    /// coupled pressure/flow consumers can factor selected response columns
    /// once and must check their accepted trajectory against their constraints.
    pub fn average_output_map(&self, output_port: usize) -> Result<(&[f64], &[f64]), PhsError> {
        let (n, m) = (self.system().n, self.system().m);
        if output_port >= m {
            return Err(PhsError::BadPortPairing);
        }
        Ok((
            &self.mean_state[output_port * n..(output_port + 1) * n],
            if self.duration == 0.0 {
                &[]
            } else {
                &self.mean_input[output_port * m..(output_port + 1) * m]
            },
        ))
    }

    /// Internal transfer view, consumed only after a successful native step.
    pub(crate) fn port_integrals(&self) -> &[Sum] { &self.scratch.port_integrals }

    /// Same-window held-input relation `y_bar[p] = free + beta*input[p]`.
    ///
    /// The selected input slot is neither read nor validated; all other held
    /// inputs and every original initial state coordinate enter `free`.
    /// This is one allocation-free O(n+m) contraction of the exact cached
    /// integral, with no trial propagation, exponential or state mutation.
    /// At zero duration it returns the original instantaneous output and zero
    /// response. A consumer requiring positive mobility must admit `beta`
    /// itself; a disconnected port can have an exactly zero response.
    pub fn average_relation(&self, initial: &[f64], input: &[f64], selected_port: usize)
        -> Result<(f64, f64), PhsError> {
        let (n, m) = (self.system().n, self.system().m);
        if selected_port >= m { return Err(PhsError::BadPortPairing); }
        if initial.len() != n || input.len() != m || !finite(initial)
            || input.iter().enumerate().any(|(p, value)| p != selected_port && !value.is_finite())
        {
            return Err(refuse("prepared ZOH mean relation requires complete finite state and other inputs"));
        }
        let state_row = &self.mean_state[selected_port*n..(selected_port+1)*n];
        let mut free = state_row.iter().zip(initial).map(|(a, x)| a*x).sum::<f64>();
        let beta = if self.duration == 0.0 { 0.0 } else {
            let input_row = &self.mean_input[selected_port*m..(selected_port+1)*m];
            free += input_row.iter().zip(input).enumerate()
                .filter(|(p, _)| *p != selected_port)
                .map(|(_, (a, u))| a*u).sum::<f64>();
            input_row[selected_port]
        };
        if !free.is_finite() {
            return Err(refuse("prepared ZOH free mean output is nonfinite"));
        }
        Ok((free, beta))
    }

    /// Evolve one held input. Heat is the original integral of grad(H)^T R grad(H),
    /// and output is its actual trajectory average, not the final port effort.
    pub fn step(&self, initial: &[f64], input: &[f64]) -> Result<StepRecord, PhsError> {
        let sys = self.system();
        let mut scratch = Scratch::new(sys.n, sys.m);
        let record = evaluate_step(sys, self.duration, self.flow.as_ref(), self.heat.as_ref(), initial, input, &mut scratch)?;
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
        self.step_with_work_inputs_into(initial, input, [], next, average_output)
            .map(|(record, _)| record)
    }

    /// Evolve `input` and independently integrate each held contribution's work
    /// `integral(work_inputs[k]^T y(t) dt)` on that same trajectory. Contributions
    /// may describe external sources and interface reactions; they do not drive
    /// additional trajectories and need not sum to `input`.
    ///
    /// Work retains the integrated port low words before outputs are rounded.
    /// Multiplying the returned mean outputs by a contribution can lose small
    /// net work between opposing ports. This call requires no warm allocation
    /// or extra propagation and leaves both output buffers unchanged on refusal.
    pub fn step_with_work_inputs_into<const K: usize>(
        &mut self, initial: &[f64], input: &[f64], work_inputs: [&[f64]; K],
        next: &mut [f64], average_output: &mut [f64],
    ) -> Result<(PreparedStepRecord, [f64; K]), PhsError> {
        let sys = self.system.get();
        if next.len() != sys.n || average_output.len() != sys.m {
            return Err(refuse("prepared ZOH output buffers must match the complete state and ports"));
        }
        let result = evaluate_step_with_work_inputs(sys, self.duration, self.flow.as_ref(),
            self.heat.as_ref(), initial, input, work_inputs, &mut self.scratch)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok(result)
    }
}

fn evaluate_step(sys: &PortHamiltonian, duration: f64, flow: Option<&PreparedQuadraticFlow>,
    heat: Option<&PositiveHeat>, initial: &[f64], input: &[f64], scratch: &mut Scratch) -> Result<PreparedStepRecord, PhsError> {
    evaluate_step_with_work_inputs(sys, duration, flow, heat, initial, input, [], scratch)
        .map(|(record, _)| record)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_step_with_work_inputs<const K: usize>(
    sys: &PortHamiltonian, duration: f64, flow: Option<&PreparedQuadraticFlow>,
    heat: Option<&PositiveHeat>, initial: &[f64], input: &[f64],
    work_inputs: [&[f64]; K], scratch: &mut Scratch,
) -> Result<(PreparedStepRecord, [f64; K]), PhsError>{
    evaluate_step_with_work_inputs_controlled(
        sys,
        duration,
        flow,
        heat,
        initial,
        input,
        work_inputs,
        scratch,
        || false,
    )
    .map_err(uncancelled)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_step_with_work_inputs_controlled<const K: usize>(
    sys: &PortHamiltonian,
    duration: f64,
    flow: Option<&PreparedQuadraticFlow>,
    heat: Option<&PositiveHeat>,
    initial: &[f64],
    input: &[f64],
    work_inputs: [&[f64]; K],
    scratch: &mut Scratch,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PreparedStepRecord, [f64; K]), PreparedStepError> {
    poll(&mut cancelled)?;
    if initial.len() != sys.n || input.len() != sys.m || !finite(initial) || !finite(input) {
        return Err(refuse("prepared ZOH requires complete finite state and input").into());
    }
    if work_inputs.iter().any(|values| values.len() != sys.m || !finite(values)) {
        return Err(refuse("prepared ZOH work observations require complete finite inputs").into());
    }
    let Some(flow) = flow else {
        scratch.effort.fill(0.0);
        sys.storage.gradient(initial, &mut scratch.effort);
        sys.output_from_effort_into(&scratch.effort, &mut scratch.output)?;
        if !finite(&scratch.output) || !sys.hamiltonian(initial).is_finite() {
            return Err(refuse("prepared ZOH initial energy/output is nonfinite").into());
        }
        scratch.endpoint.copy_from_slice(initial);
        scratch.port_integrals.fill(Sum::default());
        poll(&mut cancelled)?;
        return Ok((PreparedStepRecord { delta_h: 0.0, dissipated: 0.0, supplied: 0.0,
            newton_iters: 0, solver_residual: 0.0, solver_certificate: None }, [0.0; K]));
    };
    let (supplied, observed_work) = flow.evaluate_port_work_observed_into(initial, input,
        work_inputs, &mut scratch.augmented, &mut scratch.balanced,
        &mut scratch.endpoint, &mut scratch.integrals, Some(&mut scratch.port_integrals))?;
    poll(&mut cancelled)?;
    for (out, value) in scratch.output.iter_mut().zip(&scratch.integrals) { *out = value/duration; }
    let heat = heat.ok_or_else(||refuse("prepared positive heat is missing"))?;
    let dissipated = heat.evaluate_into(sys,initial,input,&mut scratch.effort,&mut scratch.balanced,
            &mut scratch.heat_state)?;
    poll(&mut cancelled)?;
    let delta_h = storage_energy_change(sys.storage.as_ref(), initial, &scratch.endpoint,
        &mut scratch.midpoint, &mut scratch.effort);
    let scale = energetic_scale_into(sys, initial, &mut scratch.effort)
        .max(energetic_scale_into(sys, &scratch.endpoint, &mut scratch.effort))
        .max(delta_h.abs()).max(supplied.abs()).max(dissipated.abs());
    let residual = delta_h + dissipated - supplied;
    if !finite(&scratch.output) || !sys.hamiltonian(&scratch.endpoint).is_finite()
        || !finite(&[delta_h, supplied, dissipated, scale, residual])
        || dissipated < -1e-12*scale
        || !balance_admitted(delta_h, dissipated, supplied, scale, heat.doublings()) {
        return Err(refuse("prepared ZOH physical energy balance refused").into());
    }
    poll(&mut cancelled)?;
    Ok((PreparedStepRecord { delta_h, dissipated, supplied,
        newton_iters: 0, solver_residual: 0.0, solver_certificate: None }, observed_work))
}
