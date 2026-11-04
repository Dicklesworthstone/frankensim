//! Reverse complete driven-support trajectories on a fixed accepted time grid.
//!
//! The linear motion Jacobian depends on the interval and frozen operators,
//! not on the primal state. Control recording therefore retains O(steps) times.
//! Material design optionally retains bounded full endpoint motion. Sparse
//! motion, support-reaction and stress observations seed accepted endpoints; one
//! reverse sweep returns material, force and support-control derivatives.
use super::super::material::{MaterialGradientBudget3, MaterialScaleGradient3};
use super::super::trajectory::{ControlForceKind3, ForceCotangent3};
use super::super::{TransientAdjointConfig3, TransientAdjointError3, check, number};
use super::{NodalMotionCotangent3, PrescribedTransientAdjoint3, SupportAdjointBudget3};
use crate::loads3::TetLoadCase3;
use crate::loads3::transient::events::{LoadJump3, TransientEventStep3};
use crate::loads3::transient::{TetTransientRun3, TransientStep3};
use fs_exec::Cx;

mod stress;
pub use stress::ElementStressObservation3;

/// Bounds on the accepted schedule, returned gradients and sparse observations.
#[derive(Debug, Clone, Copy)]
pub struct PrescribedRecordingBudget3 {
    /// Maximum recorded transitions, including zero-duration acceleration jumps.
    pub maximum_steps: usize,
    /// All returned q/v/a, force, support-control and thermal scalars.
    /// Times and vector metadata are instead bounded by `maximum_steps`.
    pub maximum_gradient_values: usize,
    /// Total caller-owned motion, reaction or stress rows per reverse sweep.
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
/// intervals or acceleration jumps. Left and right limits at the same physical
/// time have distinct state indices. Multiple rows are added in caller order.
/// Values have objective
/// per motion units, and are not physical displacement/velocity/acceleration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodalMotionObservation3 {
    pub state: usize,
    pub node: usize,
    pub q: [f64; 3],
    pub v: [f64; 3],
    pub a: [f64; 3],
}

/// Derivative of an objective with respect to a support-on-solid reaction [N].
/// `state == 0` is initialization; later states are accepted endpoints. `dof`
/// must be a constrained global component. Duplicate rows add in caller order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionObservation3 {
    pub state: usize,
    pub dof: usize,
    pub value: f64,
}

/// Direct derivative with respect to an external nodal force on a fixed DOF.
/// Rows retain the reaction observations' order and duplicates add. The free
/// force derivatives, including reaction-induced motion, are in `design`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionLoadCotangent3 {
    pub state: usize,
    pub dof: usize,
    pub value: f64,
}

/// Derivatives of absolute world-component prescribed displacement [objective/m].
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedDisplacementCotangent3 {
    pub time_s: f64,
    /// Sorted `RecordedPrescribedTransient3::adjoint().fixed_dofs()` order.
    pub values: Vec<f64>,
}

/// Derivatives of prescribed endpoint acceleration [objective/(m/s²)].
/// Initial prescribed acceleration remains in `PrescribedControlGradient3::initial`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedAccelerationCotangent3 {
    pub time_s: f64,
    /// IntervalEndpoint is the left endpoint of an event; RightLimit is its
    /// subsequent acceleration switch. Equal-time rows retain forward order.
    pub kind: ControlForceKind3,
    /// Sorted `RecordedPrescribedTransient3::adjoint().fixed_dofs()` order.
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedControlGradient3 {
    /// Independent initial free displacement, all initial velocities and fixed
    /// accelerations. Fixed q derivatives are in `displacements[0]`, so these
    /// entries are zero here; free a entries are zero because a0 was solved.
    pub initial: NodalMotionCotangent3,
    /// Initial load, then every recorded interval endpoint or right-limit switch.
    pub forces: Vec<ForceCotangent3>,
    /// Initial support displacement, then displacement-controlled endpoints.
    /// Omitted forward support rows still mean zero; only chain derivatives
    /// into the controls the caller actually varies.
    pub displacements: Vec<PrescribedDisplacementCotangent3>,
    /// Only acceleration-controlled endpoints and right-limit switches, in
    /// forward order. Omitted components are zero controls; chain only those
    /// the caller varies. Rows have distinct physical units from displacements.
    pub accelerations: Vec<PrescribedAccelerationCotangent3>,
    /// Derivative of the ONE frozen reduced thermal equivalent load, summed
    /// across initialization and all intervals. No material derivative implied.
    pub thermal_equivalent_force: Vec<f64>,
    pub total_adjoint_iterations: usize,
    pub maximum_adjoint_relative_residual: f64,
}

/// Opt-in storage for differentiating material parameters through a trajectory.
/// Controls alone need no motion tape; material residuals need the accepted
/// full q/v/a at initialization and each interval endpoint.
#[derive(Debug, Clone, Copy)]
pub struct PrescribedMaterialRecordingBudget3 {
    /// Maximum retained motion scalars: samples * nodes * 9. Metadata remains
    /// bounded by the recorder's independent maximum_steps allowance.
    pub maximum_motion_values: usize,
    pub gradient: MaterialGradientBudget3,
}
impl Default for PrescribedMaterialRecordingBudget3 {
    fn default() -> Self {
        Self {
            maximum_motion_values: 1_000_000,
            gradient: MaterialGradientBudget3::default(),
        }
    }
}

/// One reverse sweep for both authored controls and frozen material parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedDesignGradient3 {
    pub controls: PrescribedControlGradient3,
    /// Independent element stiffness/density scales and global Rayleigh terms.
    /// These already include stiffness-induced changes of the thermal load at
    /// fixed free thermal strain; do not chain that change a second time through
    /// controls.thermal_equivalent_force. External nodal loads are held fixed.
    pub material: MaterialScaleGradient3,
}

/// Complete fixed-grid derivatives of mixed motion and support-load objectives.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedReactionGradient3 {
    /// Includes both reaction's direct material derivative and its derivative
    /// through the original trajectory and initialization solve.
    pub design: PrescribedDesignGradient3,
    pub fixed_forces: Vec<ReactionLoadCotangent3>,
    /// Direct derivative of the one frozen full thermal equivalent force on
    /// constrained DOFs, in `adjoint().fixed_dofs()` order, summed over states.
    /// Free components are in `design.controls.thermal_equivalent_force`.
    /// Material derivatives already include stiffness-induced thermal changes;
    /// do not chain those changes through either thermal cotangent again.
    pub thermal_fixed_force: Vec<f64>,
}

enum ReverseVisit3<'a> {
    Observation {
        state: usize,
        seed: &'a mut NodalMotionCotangent3,
    },
    Residual {
        state: usize,
        force: &'a [f64],
    },
}

#[derive(Clone)]
struct MaterialTape3 {
    motion: Vec<[Vec<[f64; 3]>; 3]>,
    budget: PrescribedMaterialRecordingBudget3,
}

#[derive(Debug, Clone, Copy)]
struct RecordedSupportTransition3 {
    time_s: f64,
    h_s: f64,
    kind: SupportTransitionKind3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SupportTransitionKind3 {
    Displacement,
    Acceleration,
    AccelerationJump,
}

/// An owning recorder prevents advancement that would bypass its time record.
/// Cloning forks both the accepted original state and the recorded schedule.
#[derive(Clone)]
pub struct RecordedPrescribedTransient3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    adjoint: PrescribedTransientAdjoint3<'g, 'm, 'p, 's>,
    transitions: Vec<RecordedSupportTransition3>,
    initial_time_s: f64,
    budget: PrescribedRecordingBudget3,
    material: Option<MaterialTape3>,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Record a prescribed-support trajectory before its first finite step.
    /// The original primal and transpose policies remain independent. This
    /// derives fixed-grid motion, reaction and stress objectives. Event times,
    /// physical work reports and acceptance decisions require other pullbacks.
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
            transitions: Vec::new(),
            initial_time_s,
            budget,
            material: None,
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
        self.transitions.len() + 1
    }
    pub fn accepted_times(&self) -> impl Iterator<Item = f64> + '_ {
        std::iter::once(self.initial_time_s).chain(self.transitions.iter().map(|i| i.time_s))
    }

    /// Retain original accepted full motion for material-design derivatives.
    /// Enable this once, before advancing; storage/cancellation refusal leaves
    /// the existing recorder usable. This changes no primal numerical policy.
    pub fn enable_material_recording(
        &mut self,
        budget: PrescribedMaterialRecordingBudget3,
        cx: &Cx<'_>,
    ) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        if !self.transitions.is_empty() || self.material.is_some() {
            return Err(TransientAdjointError3::Invalid(
                "enable material recording once before the first recorded step",
            ));
        }
        admit_material_size(
            1,
            self.run.problem.nodes_m.len(),
            self.run.geometry.tets.len(),
            budget,
        )?;
        let mut motion = reserved(1)?;
        motion.push(self.run.nodal_motion(cx)?);
        check(cx)?;
        self.material = Some(MaterialTape3 { motion, budget });
        Ok(())
    }

    pub(crate) fn admit_next_step(&self, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        self.admit_transitions(1, cx)
    }

    fn admit_transitions(&self, count: usize, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        let next = self
            .transitions
            .len()
            .checked_add(count)
            .ok_or(TransientAdjointError3::Budget)?;
        admit_gradient_size(
            next,
            self.run.problem.nodes_m.len(),
            self.run.free_dofs().len(),
            self.budget,
        )?;
        if let Some(tape) = &self.material {
            let samples = next.checked_add(1).ok_or(TransientAdjointError3::Budget)?;
            admit_material_size(
                samples,
                self.run.problem.nodes_m.len(),
                self.run.geometry.tets.len(),
                tape.budget,
            )?;
        }
        Ok(())
    }

    fn reserve_transitions(
        &mut self,
        count: usize,
        cx: &Cx<'_>,
    ) -> Result<(), TransientAdjointError3> {
        self.admit_transitions(count, cx)?;
        self.transitions
            .try_reserve(count)
            .map_err(|_| TransientAdjointError3::Budget)?;
        if let Some(tape) = &mut self.material {
            tape.motion
                .try_reserve(count)
                .map_err(|_| TransientAdjointError3::Budget)?;
        }
        Ok(())
    }

    fn recorded_motion(
        &self,
        candidate: &TetTransientRun3<'_, '_, '_, '_>,
        cx: &Cx<'_>,
    ) -> Result<Option<[Vec<[f64; 3]>; 3]>, TransientAdjointError3> {
        if self.material.is_some() {
            Ok(Some(candidate.nodal_motion(cx)?))
        } else {
            Ok(None)
        }
    }

    fn append_transition(
        &mut self,
        time_s: f64,
        h_s: f64,
        kind: SupportTransitionKind3,
        motion: Option<[Vec<[f64; 3]>; 3]>,
    ) {
        self.transitions
            .push(RecordedSupportTransition3 { time_s, h_s, kind });
        if let (Some(tape), Some(motion)) = (&mut self.material, motion) {
            tape.motion.push(motion);
        }
    }

    /// Advance the original coupled solver and record its accepted time, plus
    /// full endpoint motion when material recording is enabled.
    /// Storage refusal, malformed support inputs, numerical refusal and
    /// cancellation leave both the primal and schedule reusable.
    pub fn advance_with_supports(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientAdjointError3> {
        self.advance_controlled(
            load,
            prescribed_dofs_m,
            h_s,
            SupportTransitionKind3::Displacement,
            cx,
        )
    }

    /// Integrate and record prescribed acceleration using the original coupled
    /// motion owner. Mixed displacement/acceleration intervals retain their
    /// own transpose maps and physical control units. Initial acceleration is
    /// still differentiated through the original initialization mass solve.
    pub fn advance_with_support_accelerations(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m_s2: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientAdjointError3> {
        self.advance_controlled(
            load,
            prescribed_dofs_m_s2,
            h_s,
            SupportTransitionKind3::Acceleration,
            cx,
        )
    }

    fn advance_controlled(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        values: &[(usize, f64)],
        h_s: f64,
        kind: SupportTransitionKind3,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientAdjointError3> {
        self.reserve_transitions(1, cx)?;
        let mut candidate = self.run.clone();
        let report = match kind {
            SupportTransitionKind3::Displacement => {
                candidate.advance_with_supports(load, values, h_s, cx)?
            }
            SupportTransitionKind3::Acceleration => {
                candidate.advance_with_support_accelerations(load, values, h_s, cx)?
            }
            SupportTransitionKind3::AccelerationJump => {
                return Err(TransientAdjointError3::Invalid(
                    "a finite interval requires a displacement or acceleration endpoint",
                ));
            }
        };
        let motion = self.recorded_motion(&candidate, cx)?;
        check(cx)?;
        self.append_transition(report.time_s, h_s, kind, motion);
        self.run = candidate;
        Ok(report)
    }

    /// Record a finite acceleration/load switch with continuous q/v. This
    /// consumes one recording slot and no primal time-step slot. Reverse uses
    /// the exact mass-equilibrium map; replaced acceleration has no remaining
    /// dependence on its previous value. Event times remain fixed controls.
    pub fn apply_support_acceleration_jump(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m_s2: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<LoadJump3, TransientAdjointError3> {
        self.reserve_transitions(1, cx)?;
        let mut candidate = self.run.clone();
        let report = candidate.apply_support_acceleration_jump(load, prescribed_dofs_m_s2, cx)?;
        let motion = self.recorded_motion(&candidate, cx)?;
        check(cx)?;
        self.append_transition(
            report.time_s,
            0.0,
            SupportTransitionKind3::AccelerationJump,
            motion,
        );
        self.run = candidate;
        Ok(report)
    }

    /// Atomically record a left-limit interval and its finite right-limit
    /// acceleration/load switch. Two recording slots expose both physical
    /// limits to observations. One private primal fork and optional endpoint
    /// motion samples suffice; no existing trajectory tape is copied.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_support_acceleration_jump(
        &mut self,
        load_before: &TetLoadCase3<'_, '_>,
        acceleration_before_m_s2: &[(usize, f64)],
        load_after: &TetLoadCase3<'_, '_>,
        acceleration_after_m_s2: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientEventStep3, TransientAdjointError3> {
        self.reserve_transitions(2, cx)?;
        let mut candidate = self.run.clone();
        let interval = candidate.advance_with_support_accelerations(
            load_before,
            acceleration_before_m_s2,
            h_s,
            cx,
        )?;
        let left = self.recorded_motion(&candidate, cx)?;
        let jump =
            candidate.apply_support_acceleration_jump(load_after, acceleration_after_m_s2, cx)?;
        let right = self.recorded_motion(&candidate, cx)?;
        check(cx)?;
        self.append_transition(
            interval.time_s,
            h_s,
            SupportTransitionKind3::Acceleration,
            left,
        );
        self.append_transition(
            jump.time_s,
            0.0,
            SupportTransitionKind3::AccelerationJump,
            right,
        );
        self.run = candidate;
        Ok(TransientEventStep3 { interval, jump })
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
    /// the original recurrence once, including finite acceleration switches.
    /// A weighted sampled-time integral or
    /// least-squares objective supplies its actual quadrature/residual weights
    /// here; no continuous-time or interpolation derivative is inferred.
    ///
    /// The source frame, knot times, step schedule, operators and support axes
    /// are fixed. Chain returned integrated-force and world-support rows
    /// through the caller's authored controls; direct control costs are added
    /// by their owner. Cancellation/refusal returns no partial gradient.
    pub fn reverse_with_observations(
        &self,
        terminal: &NodalMotionCotangent3,
        observations: &[NodalMotionObservation3],
        cx: &Cx<'_>,
    ) -> Result<PrescribedControlGradient3, TransientAdjointError3> {
        self.reverse_accumulating(terminal, observations, |_| Ok(()), cx)
    }

    /// Differentiate a terminal motion objective with respect to controls,
    /// per-element stiffness/density scales, and global Rayleigh damping.
    pub fn reverse_design(
        &self,
        terminal: &NodalMotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<PrescribedDesignGradient3, TransientAdjointError3> {
        self.reverse_design_with_observations(terminal, &[], cx)
    }

    /// Reverse observed motion and contract material residuals at their EXACT
    /// original accepted endpoints, including the initial mass solve. Requires
    /// enable_material_recording before the first step. Geometry, restraint
    /// axes, authored support histories, free thermal strain, and independent
    /// external nodal loads remain fixed in the material derivative.
    ///
    /// Direct objective/material or mass-dependent external-load derivatives
    /// belong to their owner. This is a fixed-grid discrete derivative, with
    /// no differentiation of solver iterations, step acceptance or event times.
    /// Failure/cancellation publishes no partial result and keeps the tape.
    pub fn reverse_design_with_observations(
        &self,
        terminal: &NodalMotionCotangent3,
        observations: &[NodalMotionObservation3],
        cx: &Cx<'_>,
    ) -> Result<PrescribedDesignGradient3, TransientAdjointError3> {
        check(cx)?;
        let tape = self
            .material
            .as_ref()
            .ok_or(TransientAdjointError3::Invalid(
                "material design needs motion recording enabled before advancement",
            ))?;
        let mut material = MaterialScaleGradient3 {
            element_scales: filled(self.run.geometry.tets.len(), [0.0; 2])?,
            rayleigh: [0.0; 2],
        };
        let controls = self.reverse_accumulating(
            terminal,
            observations,
            |visit| {
                if let ReverseVisit3::Residual { state, force } = visit {
                    let row = self.run.material_scale_pullback_nodal_at(
                        &tape.motion[state],
                        force,
                        tape.budget.gradient,
                        cx,
                    )?;
                    add_material(&mut material, row, cx)?;
                }
                Ok(())
            },
            cx,
        )?;
        check(cx)?;
        Ok(PrescribedDesignGradient3 { controls, material })
    }

    /// Reverse sparse support-reaction and motion objectives together through
    /// the original accepted trajectory. Requires material recording before
    /// advancement. Reaction means the support-on-solid constrained residual
    /// `M a + (eta M + zeta K) v + K q - f_thermal - f_external`.
    ///
    /// Supply derivatives of the objective with respect to these signed
    /// reactions. This method includes the direct dependence on material,
    /// damping and constrained loads, as well as all implicit motion changes.
    /// External loads and free thermal strain are independent of material.
    /// The time grid, geometry and restraint axes remain fixed. Refusal or
    /// cancellation returns no partial result and leaves the tape reusable.
    pub fn reverse_design_with_reactions(
        &self,
        terminal: &NodalMotionCotangent3,
        motion_observations: &[NodalMotionObservation3],
        reaction_observations: &[ReactionObservation3],
        cx: &Cx<'_>,
    ) -> Result<PrescribedReactionGradient3, TransientAdjointError3> {
        check(cx)?;
        let tape = self
            .material
            .as_ref()
            .ok_or(TransientAdjointError3::Invalid(
                "reaction design needs motion recording enabled before advancement",
            ))?;
        let count = motion_observations
            .len()
            .checked_add(reaction_observations.len())
            .ok_or(TransientAdjointError3::Budget)?;
        if count > self.budget.maximum_observations {
            return Err(TransientAdjointError3::Budget);
        }
        let fixed = self.adjoint.fixed_dofs();
        let extra = reaction_observations
            .len()
            .checked_add(fixed.len())
            .ok_or(TransientAdjointError3::Budget)?;
        let mut remaining = self.budget;
        remaining.maximum_gradient_values = remaining
            .maximum_gradient_values
            .checked_sub(extra)
            .ok_or(TransientAdjointError3::Budget)?;
        admit_gradient_size(
            self.transitions.len(),
            self.run.problem.nodes_m.len(),
            self.run.free_dofs().len(),
            remaining,
        )?;
        for row in reaction_observations {
            check(cx)?;
            if row.state >= self.sample_count() || fixed.binary_search(&row.dof).is_err() {
                return Err(TransientAdjointError3::Invalid(
                    "reaction observation must name an accepted state and constrained global DOF",
                ));
            }
            number(row.value)?;
        }
        let mut heads = filled(self.sample_count(), usize::MAX)?;
        let mut next = filled(reaction_observations.len(), usize::MAX)?;
        for (i, row) in reaction_observations.iter().enumerate().rev() {
            check(cx)?;
            next[i] = heads[row.state];
            heads[row.state] = i;
        }
        let mut fixed_forces = reserved(reaction_observations.len())?;
        let mut thermal_fixed_force = filled(fixed.len(), 0.0)?;
        for row in reaction_observations {
            check(cx)?;
            // Validated above; retain caller row order without hiding duplicate
            // terms or pretending constrained forces drive free equilibrium.
            let index = fixed
                .binary_search(&row.dof)
                .map_err(|_| TransientAdjointError3::Invalid("constrained reaction DOF"))?;
            fixed_forces.push(ReactionLoadCotangent3 {
                state: row.state,
                dof: row.dof,
                value: -row.value,
            });
            thermal_fixed_force[index] = number(thermal_fixed_force[index] - row.value)?;
        }
        let mut material = MaterialScaleGradient3 {
            element_scales: filled(self.run.geometry.tets.len(), [0.0; 2])?,
            rayleigh: [0.0; 2],
        };
        let mut weights = filled(fixed.len(), 0.0)?;
        let mut rows = reserved(fixed.len())?;
        let controls = self.reverse_accumulating(
            terminal,
            motion_observations,
            |visit| {
                match visit {
                    ReverseVisit3::Observation { state, seed } => {
                        if heads[state] == usize::MAX {
                            return Ok(());
                        }
                        weights.fill(0.0);
                        rows.clear();
                        let mut index = heads[state];
                        while index != usize::MAX {
                            check(cx)?;
                            let row = reaction_observations[index];
                            let slot = fixed.binary_search(&row.dof).map_err(|_| {
                                TransientAdjointError3::Invalid("constrained reaction DOF")
                            })?;
                            weights[slot] = number(weights[slot] + row.value)?;
                            index = next[index];
                        }
                        for (&dof, &value) in fixed.iter().zip(&weights) {
                            check(cx)?;
                            if value != 0.0 {
                                rows.push((dof, value));
                            }
                        }
                        if rows.is_empty() {
                            return Ok(());
                        }
                        let motion = self.adjoint.reaction_motion_pullback(&rows, cx)?;
                        for field in [&mut seed.q, &mut seed.v, &mut seed.a]
                            .into_iter()
                            .zip([motion.q, motion.v, motion.a])
                        {
                            for (sum, value) in field.0.iter_mut().zip(field.1) {
                                check(cx)?;
                                for c in 0..3 {
                                    sum[c] = number(sum[c] + value[c])?;
                                }
                            }
                        }
                        let direct = self.run.reaction_material_pullback_nodal_at(
                            &tape.motion[state],
                            &rows,
                            tape.budget.gradient,
                            cx,
                        )?;
                        add_material(&mut material, direct, cx)?;
                    }
                    ReverseVisit3::Residual { state, force } => {
                        let implicit = self.run.material_scale_pullback_nodal_at(
                            &tape.motion[state],
                            force,
                            tape.budget.gradient,
                            cx,
                        )?;
                        add_material(&mut material, implicit, cx)?;
                    }
                }
                Ok(())
            },
            cx,
        )?;
        check(cx)?;
        Ok(PrescribedReactionGradient3 {
            design: PrescribedDesignGradient3 { controls, material },
            fixed_forces,
            thermal_fixed_force,
        })
    }

    fn reverse_accumulating(
        &self,
        terminal: &NodalMotionCotangent3,
        observations: &[NodalMotionObservation3],
        mut visit: impl FnMut(ReverseVisit3<'_>) -> Result<(), TransientAdjointError3>,
        cx: &Cx<'_>,
    ) -> Result<PrescribedControlGradient3, TransientAdjointError3> {
        check(cx)?;
        let nodes = self.run.problem.nodes_m.len();
        let free = self.run.free_dofs().len();
        terminal.validate(nodes, cx)?;
        admit_gradient_size(self.transitions.len(), nodes, free, self.budget)?;
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
        let acceleration_count = self
            .transitions
            .iter()
            .filter(|step| step.kind != SupportTransitionKind3::Displacement)
            .count();
        let mut displacements = reserved(self.sample_count() - acceleration_count)?;
        let mut accelerations = reserved(acceleration_count)?;
        let mut thermal = filled(free, 0.0)?;
        let mut iterations = 0usize;
        let mut maximum_relative = 0.0_f64;
        inject(
            &mut seed,
            self.transitions.len(),
            observations,
            &heads,
            &next,
            cx,
        )?;
        visit(ReverseVisit3::Observation {
            state: self.transitions.len(),
            seed: &mut seed,
        })?;
        for (i, interval) in self.transitions.iter().enumerate().rev() {
            check(cx)?;
            let (previous, force, control, solve_iterations, relative) = match interval.kind {
                SupportTransitionKind3::Displacement => {
                    let pb = self.adjoint.pullback_step(interval.h_s, &seed, cx)?;
                    (
                        pb.previous,
                        pb.force,
                        pb.prescribed_displacement,
                        pb.iterations,
                        pb.true_relative_residual,
                    )
                }
                SupportTransitionKind3::Acceleration | SupportTransitionKind3::AccelerationJump => {
                    let pb = if interval.kind == SupportTransitionKind3::AccelerationJump {
                        self.adjoint.pullback_acceleration_jump(&seed, cx)?
                    } else {
                        self.adjoint
                            .pullback_acceleration_step(interval.h_s, &seed, cx)?
                    };
                    (
                        pb.previous,
                        pb.force,
                        pb.prescribed_acceleration,
                        pb.iterations,
                        pb.true_relative_residual,
                    )
                }
            };
            let kind = if interval.kind == SupportTransitionKind3::AccelerationJump {
                ControlForceKind3::RightLimit
            } else {
                ControlForceKind3::IntervalEndpoint
            };
            visit(ReverseVisit3::Residual {
                state: i + 1,
                force: &force,
            })?;
            add_thermal(&mut thermal, &force, cx)?;
            iterations = iterations
                .checked_add(solve_iterations)
                .ok_or(TransientAdjointError3::Budget)?;
            maximum_relative = maximum_relative.max(relative);
            forces.push(ForceCotangent3 {
                time_s: interval.time_s,
                kind,
                values: force,
            });
            if interval.kind == SupportTransitionKind3::Displacement {
                displacements.push(PrescribedDisplacementCotangent3 {
                    time_s: interval.time_s,
                    values: control,
                });
            } else {
                accelerations.push(PrescribedAccelerationCotangent3 {
                    time_s: interval.time_s,
                    kind,
                    values: control,
                });
            }
            seed = previous;
            inject(&mut seed, i, observations, &heads, &next, cx)?;
            visit(ReverseVisit3::Observation {
                state: i,
                seed: &mut seed,
            })?;
        }
        let mut initial = self.adjoint.pullback_initialization(&seed, cx)?;
        visit(ReverseVisit3::Residual {
            state: 0,
            force: &initial.force,
        })?;
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
        accelerations.reverse();
        check(cx)?;
        Ok(PrescribedControlGradient3 {
            initial: initial.initial,
            forces,
            displacements,
            accelerations,
            thermal_equivalent_force: thermal,
            total_adjoint_iterations: iterations,
            maximum_adjoint_relative_residual: maximum_relative,
        })
    }
}

fn admit_material_size(
    samples: usize,
    nodes: usize,
    elements: usize,
    budget: PrescribedMaterialRecordingBudget3,
) -> Result<(), TransientAdjointError3> {
    let full = nodes.checked_mul(3).ok_or(TransientAdjointError3::Budget)?;
    let values = full
        .checked_mul(3)
        .and_then(|n| n.checked_mul(samples))
        .ok_or(TransientAdjointError3::Budget)?;
    if values > budget.maximum_motion_values
        || full > budget.gradient.maximum_full_dofs
        || elements > budget.gradient.maximum_elements
    {
        return Err(TransientAdjointError3::Budget);
    }
    Ok(())
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
fn add_material(
    target: &mut MaterialScaleGradient3,
    source: MaterialScaleGradient3,
    cx: &Cx<'_>,
) -> Result<(), TransientAdjointError3> {
    for (sum, value) in target.element_scales.iter_mut().zip(source.element_scales) {
        check(cx)?;
        for c in 0..2 {
            sum[c] = number(sum[c] + value[c])?;
        }
    }
    for c in 0..2 {
        target.rayleigh[c] = number(target.rayleigh[c] + source.rayleigh[c])?;
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
