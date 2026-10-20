//! Fixed-grid scan adjoints with short replay tapes and caloric checkpoints.
//!
//! A forward pass stores only segment-boundary states. Reverse processing restores
//! each segment, rebuilds its original implicit tangents, and drops them before
//! the preceding segment. This costs one additional primal evaluation per
//! endpoint on success; it is uniform segmentation, not optimal Revolve.
//! Restarts stay private to this consumed trial run, never a public operation
//! for replacing live physical history or refunding accepted application work.

use super::{ScanObjectiveBudget, ScanObjectiveEvaluation, ScanObservation,
    ScanObjectiveSampling, ScanSensitivityBudget, ScanSensitivityError,
    ScanStepLinearization, ScanningEnthalpyRun, ScanStep, EnthalpyStepConfig,
    LinearConfig, ReverseState, admit_layout, advance_binding, observed, finite,
    invalid, poll, zeros};
use fs_exec::Cx;

mod responses;
pub use responses::{CheckpointedScanResponses, evaluate_scan_objective_checkpointed_with_peaks};
pub(crate) use responses::evaluate_observations_checkpointed;

/// Separate checkpoint storage and recomputation allowances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanCheckpointBudget {
    /// At most this many endpoint linearizations coexist. Must be nonzero.
    pub segment_steps: usize,
    /// (ceil(endpoints / segment_steps) + 1) times mesh vertices. Each checkpoint
    /// holds THREE arrays (h, T, liquid fraction), a time and a logical step count.
    /// This is a shape limit, not bytes, allocator overhead or a total memory cap.
    pub max_checkpoint_vertices: usize,
    /// Original plus recomputed primal endpoints; success needs 2*endpoints.
    /// Inner Newton/Krylov/assembly budgets remain separate and unchanged.
    pub max_forward_solves: usize,
}

/// Actual successful work and retained shapes, not an estimated speedup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanCheckpointWork {
    pub original_forward_solves: usize,
    pub recomputed_forward_solves: usize,
    pub checkpoint_states: usize,
    pub checkpoint_vertices: usize,
    pub peak_retained_linearizations: usize,
}

/// Original physical objective/adjoint, plus its memory/work tradeoff.
#[derive(Debug)]
pub struct CheckpointedScanObjective {
    pub objective: ScanObjectiveEvaluation,
    pub work: ScanCheckpointWork,
}

/// Evaluate a scan objective while retaining only a short sparse-tangent tape.
///
/// `budget.max_steps` still bounds the WHOLE trajectory. Its retained vertex,
/// source-vertex and radiation-entry allowances now bound ONE replay segment.
/// Checkpoint arrays have their own allowance. No solver tolerance, mesh,
/// event side or time grid changes. Segment endpoints and fitting observations
/// must replay bit-exactly before their derivatives are used. This is a runtime
/// determinism check, not an interval certificate or durable checkpoint format.
///
/// # Errors
/// Original source/phase/material/solver/energy/transpose refusals, invalid
/// budgets, omitted events/motion stops, replay mismatch or cancellation. The
/// trial run is consumed; neither a partial gradient nor a changed caller-owned
/// trajectory escapes. Logical counts are restored only in this private trial;
/// successful original/recomputed solves are reported separately.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_objective_checkpointed(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    checkpoints: ScanCheckpointBudget,
) -> Result<CheckpointedScanObjective, ScanSensitivityError> {
    let result = evaluate_scan_objective_checkpointed_with_peaks(cx, run, endpoints,
        observations, sampling, step, transpose, budget, &[], &[],
        super::responses::ScanResponseBudget { max_responses: 0, max_gradient_entries: 0,
            max_transpose_solves: endpoints.len() }, 0, checkpoints)?;
    Ok(CheckpointedScanObjective { objective: result.objective, work: result.work })
}

fn admit_checkpoints(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    count: usize, budget: ScanCheckpointBudget,
) -> Result<usize, ScanSensitivityError> {
    poll(cx, 0)?;
    if count == 0 || budget.segment_steps == 0 {
        return Err(invalid("checkpointed scan needs nonempty endpoints and a positive segment length"));
    }
    let segment = budget.segment_steps.min(count);
    let states = count.div_ceil(segment).checked_add(1)
        .ok_or_else(|| invalid("scan checkpoint count overflow"))?;
    if states.checked_mul(run.conduction().temperature().len())
        .is_none_or(|n| n > budget.max_checkpoint_vertices)
    { return Err(invalid("scan caloric checkpoints exceed their separate retained vertex allowance")); }
    if count.checked_mul(2).is_none_or(|n| n > budget.max_forward_solves) {
        return Err(invalid("scan original and recomputed endpoints exceed their forward-solve allowance"));
    }
    Ok(segment)
}

// Complete evolving state. Geometry, charts, reference masses, radiation,
// contacts, sources and motion policy stay in the SAME immutable run.
struct Snapshot {
    h: Vec<f64>,
    temperature: Vec<f64>,
    fraction: Vec<f64>,
    time_s: f64,
    accepted_steps: usize,
}
impl Snapshot {
    fn capture(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>)
        -> Result<Self, ScanSensitivityError>
    {
        let state = run.conduction();
        Ok(Self { h: copy(cx, state.specific_enthalpy_j_kg())?,
            temperature: copy(cx, state.temperature())?,
            fraction: copy(cx, state.liquid_mass_fraction())?,
            time_s: state.time_s(), accepted_steps: state.accepted_steps() })
    }

    fn restore(&self, cx: &Cx<'_>, run: &mut ScanningEnthalpyRun<'_, '_, '_>)
        -> Result<(), ScanSensitivityError>
    {
        let h = copy(cx, &self.h)?;
        let temperature = copy(cx, &self.temperature)?;
        let fraction = copy(cx, &self.fraction)?;
        poll(cx, self.accepted_steps)?;
        run.run.specific_h = h;
        run.run.temperature = temperature;
        run.run.liquid_mass_fraction = fraction;
        run.run.time_s = self.time_s;
        run.run.accepted_steps = self.accepted_steps;
        Ok(())
    }

    fn check(&self, cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>)
        -> Result<(), ScanSensitivityError>
    {
        let state = run.conduction();
        same_scalar(self.time_s, state.time_s())?;
        if self.accepted_steps != state.accepted_steps() {
            return Err(invalid("scan checkpoint replay changed its logical step count"));
        }
        for (expected, actual) in [(&self.h[..], state.specific_enthalpy_j_kg()),
            (&self.temperature[..], state.temperature()), (&self.fraction[..], state.liquid_mass_fraction())]
        {
            if expected.len() != actual.len() { return Err(invalid("scan checkpoint state shape changed")); }
            for (i, (&a, &b)) in expected.iter().zip(actual).enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                same_scalar(a, b)?;
            }
        }
        poll(cx, self.accepted_steps)
    }
}

// Solver-specific replay owner, not a second numerical integrator or a public
// rewind operation. Snapshot allocation is admitted before the forward pass.
struct ReplayedScan<'h, 'm, 'c> {
    run: ScanningEnthalpyRun<'h, 'm, 'c>,
    snapshots: Vec<Snapshot>,
    segment_steps: usize,
    final_step: ScanStep,
    work: ScanCheckpointWork,
}
impl<'h, 'm, 'c> ReplayedScan<'h, 'm, 'c> {
    #[allow(clippy::too_many_arguments)]
    fn forward(cx: &Cx<'_>, mut run: ScanningEnthalpyRun<'h, 'm, 'c>, endpoints: &[f64],
        sampling: ScanObjectiveSampling, step: EnthalpyStepConfig, segment_steps: usize,
        mut observe: impl FnMut(usize, &ScanStep) -> Result<(), ScanSensitivityError>,
    ) -> Result<Self, ScanSensitivityError> {
        let states = endpoints.len().div_ceil(segment_steps) + 1;
        let mut snapshots = Vec::new();
        snapshots.try_reserve_exact(states).map_err(|_| invalid("scan checkpoint allocation refused"))?;
        snapshots.push(Snapshot::capture(cx, &run)?);
        let mut final_step = None;
        for (index, &time) in endpoints.iter().enumerate() {
            poll(cx, index)?;
            let report = match sampling {
                ScanObjectiveSampling::RightEndpoint => run.advance_to(cx, time, step)?,
                ScanObjectiveSampling::EventLeft => {
                    let result = run.advance_event_step(cx, time, time-run.conduction().time_s(), step)?;
                    if result.end_time_s != time {
                        return Err(invalid("checkpointed scan grid omits a source/boundary event or selected motion stop"));
                    }
                    result.report
                }
            };
            observe(index, &report)?;
            if (index+1) % segment_steps == 0 || index+1 == endpoints.len() {
                snapshots.push(Snapshot::capture(cx, &run)?);
            }
            if index+1 == endpoints.len() { final_step = Some(report); }
        }
        let vertices = run.conduction().temperature().len();
        poll(cx, endpoints.len())?;
        Ok(Self { run, snapshots, segment_steps,
            final_step: final_step.ok_or_else(|| invalid("missing checkpointed scan final report"))?,
            work: ScanCheckpointWork { original_forward_solves: endpoints.len(),
                recomputed_forward_solves: 0, checkpoint_states: states,
                checkpoint_vertices: states*vertices, peak_retained_linearizations: 0 } })
    }

    fn segments(&self) -> usize { self.snapshots.len()-1 }

    fn segment(&mut self, cx: &Cx<'_>, index: usize, endpoints: &[f64],
        sampling: ScanObjectiveSampling, step: EnthalpyStepConfig, per_step: ScanSensitivityBudget,
    ) -> Result<(usize, Vec<ScanStepLinearization<'m>>), ScanSensitivityError> {
        let begin = index*self.segment_steps;
        let end = begin + self.segment_steps.min(endpoints.len()-begin);
        self.snapshots[index].restore(cx, &mut self.run)?;
        let mut tape = Vec::new();
        tape.try_reserve_exact(end-begin).map_err(|_| invalid("scan replay tape allocation refused"))?;
        for &time in &endpoints[begin..end] {
            poll(cx, self.work.recomputed_forward_solves)?;
            let binding = advance_binding(cx, &mut self.run, time, sampling, step, per_step)?;
            self.work.recomputed_forward_solves += 1;
            tape.push(binding);
        }
        self.snapshots[index+1].check(cx, &self.run)?;
        self.work.peak_retained_linearizations = self.work.peak_retained_linearizations.max(tape.len());
        Ok((begin, tape))
    }
}

fn copy(cx: &Cx<'_>, input: &[f64]) -> Result<Vec<f64>, ScanSensitivityError> {
    let mut output = zeros(input.len())?;
    for (i, (out, &value)) in output.iter_mut().zip(input).enumerate() {
        if i % 256 == 0 { poll(cx, i)?; }
        *out = value;
    }
    Ok(output)
}

fn same_scalar(expected: f64, actual: f64) -> Result<(), ScanSensitivityError> {
    if expected.to_bits() != actual.to_bits() {
        return Err(invalid("scan checkpoint replay changed an accepted field or observation; no gradient is published"));
    }
    Ok(())
}
