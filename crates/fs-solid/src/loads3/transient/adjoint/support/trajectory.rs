//! Reverse complete driven-support trajectories on a fixed accepted time grid.
//!
//! The linear motion Jacobian depends on the interval and frozen operators,
//! not on the primal state. Recording therefore retains O(steps) times rather
//! than a full motion tape. Sparse observations seed any accepted endpoint;
//! one reverse sweep returns every force and support-displacement derivative.
use super::super::trajectory::{ControlForceKind3, ForceCotangent3};
use super::super::{TransientAdjointConfig3, TransientAdjointError3, check, number};
use super::{NodalMotionCotangent3, PrescribedTransientAdjoint3, SupportAdjointBudget3};
use crate::loads3::TetLoadCase3;
use crate::loads3::transient::{TetTransientRun3, TransientStep3};
use fs_exec::Cx;

/// Bounds on the accepted schedule, returned gradients and sparse observations.
#[derive(Debug, Clone, Copy)]
pub struct PrescribedRecordingBudget3 {
    pub maximum_steps: usize,
    /// All returned q/v/a, force, displacement-control and thermal scalars.
    /// Times and vector metadata are instead bounded by `maximum_steps`.
    pub maximum_gradient_values: usize,
    /// Number of caller-owned node/endpoint observation rows per reverse sweep.
    pub maximum_observations: usize,
}
impl Default for PrescribedRecordingBudget3 {
    fn default() -> Self {
        Self {
            maximum_steps: 10_000,
            maximum_gradient_values: 1_000_000,
            maximum_observations: 100_000,
        }
    }
}

/// A scalar objective's derivative at one node of an accepted motion state.
/// `state == 0` addresses initialization; subsequent states follow accepted
/// intervals. Multiple rows are added in caller order. Values have objective
/// per motion units, and are not physical displacement/velocity/acceleration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodalMotionObservation3 {
    pub state: usize,
    pub node: usize,
    pub q: [f64; 3],
    pub v: [f64; 3],
    pub a: [f64; 3],
}

/// Derivatives of absolute world-component prescribed displacement [objective/m].
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedDisplacementCotangent3 {
    pub time_s: f64,
    /// Sorted `RecordedPrescribedTransient3::adjoint().fixed_dofs()` order.
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedControlGradient3 {
    /// Independent initial free displacement, all initial velocities and fixed
    /// accelerations. Fixed q derivatives are in `displacements[0]`, so these
    /// entries are zero here; free a entries are zero because a0 was solved.
    pub initial: NodalMotionCotangent3,
    /// Initial load, then accepted interval endpoints in forward time order.
    pub forces: Vec<ForceCotangent3>,
    /// Initial support displacement, then every accepted interval endpoint.
    /// Omitted forward support rows still mean zero; only chain derivatives
    /// into the controls the caller actually varies.
    pub displacements: Vec<PrescribedDisplacementCotangent3>,
    /// Derivative of the ONE frozen reduced thermal equivalent load, summed
    /// across initialization and all intervals. No material derivative implied.
    pub thermal_equivalent_force: Vec<f64>,
    pub total_adjoint_iterations: usize,
    pub maximum_adjoint_relative_residual: f64,
}

#[derive(Debug, Clone, Copy)]
struct Interval3 {
    time_s: f64,
    h_s: f64,
}

/// An owning recorder prevents advancement that would bypass its time record.
/// Cloning forks both the accepted original state and the recorded schedule.
#[derive(Clone)]
pub struct RecordedPrescribedTransient3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    adjoint: PrescribedTransientAdjoint3<'g, 'm, 'p, 's>,
    intervals: Vec<Interval3>,
    initial_time_s: f64,
    budget: PrescribedRecordingBudget3,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Record a prescribed-support trajectory before its first finite step.
    /// The original primal and transpose policies remain independent. This
    /// derives fixed-grid motion objectives, not event times, reactions,
    /// stresses, physical work reports or derivatives of acceptance decisions.
    pub fn record_prescribed_controls(
        self,
        budget: PrescribedRecordingBudget3,
        adjoint_config: TransientAdjointConfig3,
        support_budget: SupportAdjointBudget3,
        cx: &Cx<'_>,
    ) -> Result<RecordedPrescribedTransient3<'g, 'm, 'p, 's>, TransientAdjointError3> {
        check(cx)?;
        if self.state().step != 0 {
            return Err(TransientAdjointError3::Invalid(
                "start prescribed control recording before the first finite step",
            ));
        }
        admit_gradient_size(
            0,
            self.problem.nodes_m.len(),
            self.free_dofs().len(),
            budget,
        )?;
        let adjoint = self.prescribed_adjoint(adjoint_config, support_budget, cx)?;
        let initial_time_s = self.state().time_s;
        check(cx)?;
        Ok(RecordedPrescribedTransient3 {
            run: self,
            adjoint,
            intervals: Vec::new(),
            initial_time_s,
            budget,
        })
    }
}

impl<'g, 'm, 'p, 's> RecordedPrescribedTransient3<'g, 'm, 'p, 's> {
    pub fn run(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> {
        &self.run
    }
    pub fn adjoint(&self) -> &PrescribedTransientAdjoint3<'g, 'm, 'p, 's> {
        &self.adjoint
    }
    pub fn sample_count(&self) -> usize {
        self.intervals.len() + 1
    }
    pub fn accepted_times(&self) -> impl Iterator<Item = f64> + '_ {
        std::iter::once(self.initial_time_s).chain(self.intervals.iter().map(|i| i.time_s))
    }

    pub(crate) fn admit_next_step(&self, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        let next = self
            .intervals
            .len()
            .checked_add(1)
            .ok_or(TransientAdjointError3::Budget)?;
        admit_gradient_size(
            next,
            self.run.problem.nodes_m.len(),
            self.run.free_dofs().len(),
            self.budget,
        )
    }

    /// Advance the original coupled solver and record only its accepted time.
    /// Storage refusal, malformed support inputs, numerical refusal and
    /// cancellation leave both the primal and schedule reusable.
    pub fn advance_with_supports(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientAdjointError3> {
        self.admit_next_step(cx)?;
        self.intervals
            .try_reserve(1)
            .map_err(|_| TransientAdjointError3::Budget)?;
        let mut candidate = self.run.clone();
        let report = candidate.advance_with_supports(load, prescribed_dofs_m, h_s, cx)?;
        check(cx)?;
        self.intervals.push(Interval3 {
            time_s: report.time_s,
            h_s,
        });
        self.run = candidate;
        Ok(report)
    }

    /// Reverse a terminal full-nodal motion objective through all intervals.
    pub fn reverse(
        &self,
        terminal: &NodalMotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<PrescribedControlGradient3, TransientAdjointError3> {
        self.reverse_with_observations(terminal, &[], cx)
    }

    /// Add sparse objective seeds at arbitrary accepted endpoints, then reverse
    /// the original recurrence once. A weighted sampled-time integral or
    /// least-squares objective supplies its actual quadrature/residual weights
    /// here; no continuous-time or interpolation derivative is inferred.
    ///
    /// The source frame, knot times, step schedule, operators and support axes
    /// are fixed. Chain returned integrated-force and world-displacement rows
    /// through the caller's authored controls; direct control costs are added
    /// by their owner. Cancellation/refusal returns no partial gradient.
    pub fn reverse_with_observations(
        &self,
        terminal: &NodalMotionCotangent3,
        observations: &[NodalMotionObservation3],
        cx: &Cx<'_>,
    ) -> Result<PrescribedControlGradient3, TransientAdjointError3> {
        check(cx)?;
        let nodes = self.run.problem.nodes_m.len();
        let free = self.run.free_dofs().len();
        terminal.validate(nodes, cx)?;
        admit_gradient_size(self.intervals.len(), nodes, free, self.budget)?;
        if observations.len() > self.budget.maximum_observations {
            return Err(TransientAdjointError3::Budget);
        }
        for observation in observations {
            check(cx)?;
            if observation.state >= self.sample_count() || observation.node >= nodes {
                return Err(TransientAdjointError3::Invalid(
                    "observation must name an accepted state and existing node",
                ));
            }
            for value in observation
                .q
                .into_iter()
                .chain(observation.v)
                .chain(observation.a)
            {
                number(value)?;
            }
        }
        // Linked indices give O(steps + observations) work and preserve caller
        // summation order without storing a full-nodal seed at every endpoint.
        let mut heads = filled(self.sample_count(), usize::MAX)?;
        let mut next = filled(observations.len(), usize::MAX)?;
        for (i, observation) in observations.iter().enumerate().rev() {
            check(cx)?;
            next[i] = heads[observation.state];
            heads[observation.state] = i;
        }
        let mut seed = terminal.clone();
        let mut forces = reserved(self.sample_count())?;
        let mut displacements = reserved(self.sample_count())?;
        let mut thermal = filled(free, 0.0)?;
        let mut iterations = 0usize;
        let mut maximum_relative = 0.0_f64;
        inject(
            &mut seed,
            self.intervals.len(),
            observations,
            &heads,
            &next,
            cx,
        )?;
        for (i, interval) in self.intervals.iter().enumerate().rev() {
            check(cx)?;
            let pb = self.adjoint.pullback_step(interval.h_s, &seed, cx)?;
            add_thermal(&mut thermal, &pb.force, cx)?;
            iterations = iterations
                .checked_add(pb.iterations)
                .ok_or(TransientAdjointError3::Budget)?;
            maximum_relative = maximum_relative.max(pb.true_relative_residual);
            forces.push(ForceCotangent3 {
                time_s: interval.time_s,
                kind: ControlForceKind3::IntervalEndpoint,
                values: pb.force,
            });
            displacements.push(PrescribedDisplacementCotangent3 {
                time_s: interval.time_s,
                values: pb.prescribed_displacement,
            });
            seed = pb.previous;
            inject(&mut seed, i, observations, &heads, &next, cx)?;
        }
        let mut initial = self.adjoint.pullback_initialization(&seed, cx)?;
        add_thermal(&mut thermal, &initial.force, cx)?;
        iterations = iterations
            .checked_add(initial.iterations)
            .ok_or(TransientAdjointError3::Budget)?;
        maximum_relative = maximum_relative.max(initial.true_relative_residual);
        let mut initial_displacement = reserved(self.adjoint.fixed_dofs().len())?;
        for &dof in self.adjoint.fixed_dofs() {
            check(cx)?;
            initial_displacement.push(initial.initial.q[dof / 3][dof % 3]);
            initial.initial.q[dof / 3][dof % 3] = 0.0;
        }
        forces.push(ForceCotangent3 {
            time_s: self.initial_time_s,
            kind: ControlForceKind3::Initial,
            values: initial.force,
        });
        displacements.push(PrescribedDisplacementCotangent3 {
            time_s: self.initial_time_s,
            values: initial_displacement,
        });
        forces.reverse();
        displacements.reverse();
        check(cx)?;
        Ok(PrescribedControlGradient3 {
            initial: initial.initial,
            forces,
            displacements,
            thermal_equivalent_force: thermal,
            total_adjoint_iterations: iterations,
            maximum_adjoint_relative_residual: maximum_relative,
        })
    }
}

fn admit_gradient_size(
    steps: usize,
    nodes: usize,
    free: usize,
    budget: PrescribedRecordingBudget3,
) -> Result<(), TransientAdjointError3> {
    let full = nodes.checked_mul(3).ok_or(TransientAdjointError3::Budget)?;
    let values = steps
        .checked_add(1)
        .and_then(|n| n.checked_mul(full))
        .and_then(|n| {
            full.checked_mul(3)
                .and_then(|initial| n.checked_add(initial))
        })
        .and_then(|n| n.checked_add(free))
        .ok_or(TransientAdjointError3::Budget)?;
    if steps > budget.maximum_steps || values > budget.maximum_gradient_values {
        return Err(TransientAdjointError3::Budget);
    }
    Ok(())
}

fn reserved<T>(capacity: usize) -> Result<Vec<T>, TransientAdjointError3> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| TransientAdjointError3::Budget)?;
    Ok(values)
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, TransientAdjointError3> {
    let mut values = reserved(count)?;
    values.resize(count, value);
    Ok(values)
}
fn add_thermal(
    target: &mut [f64],
    source: &[f64],
    cx: &Cx<'_>,
) -> Result<(), TransientAdjointError3> {
    for (i, (sum, &value)) in target.iter_mut().zip(source).enumerate() {
        if i % 256 == 0 {
            check(cx)?;
        }
        *sum = number(*sum + value)?;
    }
    Ok(())
}
fn inject(
    seed: &mut NodalMotionCotangent3,
    state: usize,
    observations: &[NodalMotionObservation3],
    heads: &[usize],
    next: &[usize],
    cx: &Cx<'_>,
) -> Result<(), TransientAdjointError3> {
    let mut index = heads[state];
    while index != usize::MAX {
        check(cx)?;
        let row = observations[index];
        for c in 0..3 {
            seed.q[row.node][c] = number(seed.q[row.node][c] + row.q[c])?;
            seed.v[row.node][c] = number(seed.v[row.node][c] + row.v[c])?;
            seed.a[row.node][c] = number(seed.a[row.node][c] + row.a[c])?;
        }
        index = next[index];
    }
    Ok(())
}
