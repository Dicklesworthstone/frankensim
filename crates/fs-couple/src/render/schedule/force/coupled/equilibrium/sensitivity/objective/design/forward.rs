//! Forward observations of the same immutable multi-case mechanical problem.
//! No tangent or adjoint is built. Numerical preload checks still apply, but
//! a physical displacement does not require differentiability at contact onset.
use super::*;

mod responses;
pub use responses::{ForwardConstraintResult, ForwardDesignResponses};

// Shared preparation belongs to the physical evaluator, not its consumers.
pub(super) struct PreparedDesign {
    pub parameters: Vec<f64>,
    pub springs: Vec<ModalConnection>,
    pub contacts: Vec<(ModalContact, ModalContactConfig)>,
    pub cases: Vec<DesignLoadCase>,
}
pub(super) struct SolvedDesignCase {
    pub network: CoupledModalSystem,
    pub external: Vec<f64>,
    pub actuators: Vec<ModalAttachment>,
    pub stored_energy_j: f64,
}

/// Complete primal observations for one independent experiment.
/// This type deliberately has no gradient or derivative-admission report.
#[derive(Clone, Debug, PartialEq)]
pub struct ForwardDesignCase {
    /// Physical displacement [m] at every declared target attachment, in order.
    /// The objective's targets, scales and weights do not change observations.
    pub observations_m: Vec<f64>,
    /// Stored component, bilateral spring and contact energy [J] returned by
    /// the original preload solver. Not an inferred actuator-work history.
    pub stored_energy_j: f64,
}

/// Every independent case, or no result on error. Absence of derivatives is
/// explicit in the type, not a fabricated zero gradient or adjoint residual.
#[derive(Clone, Debug, PartialEq)]
pub struct ForwardDesignEvaluation {
    /// Decoded physical parameter values, in variable declaration order.
    pub physical_parameters: Vec<f64>,
    /// Complete observations in load-case order.
    pub cases: Vec<ForwardDesignCase>,
    /// Number of declared response constraints NOT evaluated by this query.
    /// Observations alone must not be interpreted as constraint feasibility.
    pub unassessed_response_constraints: usize,
}

impl EquilibriumDesign {
    /// Evaluate physical observations without building tangents or adjoints.
    ///
    /// Uses exactly the same candidate binding, force projection, independent
    /// zero-state templates, preload algorithms and cumulative work accounting
    /// as `evaluate`. The preload's setup, force, penetration, full-network
    /// residual and energy checks are unchanged. An original physical failure
    /// is terminal for this call; no failed case is replaced or dropped.
    ///
    /// The derivative-only `SensitivityBudget` does not restrict this query:
    /// observation work is bounded by the admitted modal/port/case counts.
    /// Contact onset and points inside a derivative exclusion margin may be
    /// observed when the original primal solver succeeds. This grants NO
    /// derivative, stability, interval, physical-validation or global-uniqueness
    /// claim. Unsupported contact-only restraint of free modes still refuses.
    ///
    /// Cancellation is checked before cases and observations and in the same
    /// preload owners. A failure never refunds spent work or mutates templates.
    pub fn evaluate_forward(&self, point: &[f64], control: &mut DesignControl, gate: &CancelGate)
        -> Result<ForwardDesignEvaluation, DesignError>
    {
        self.forward_query(point, control, gate, false).map(|result| result.observations)
    }

    /// Evaluate EVERY authored response constraint on the same solved cases as
    /// the observations, without constructing an adjoint or requiring an activity
    /// margin. Constraint values are signed physical observations; a violated
    /// design requirement is not a failed physical solve. Equality residuals
    /// carry no implicit acceptance tolerance. No derivative or certificate is
    /// returned, and no partially evaluated family escapes on any refusal.
    ///
    /// Reuses the original spring reaction and fs-dcontact potential gradient.
    /// Work is bounded by admitted modes/cases/ports and the existing 64-row
    /// constraint cap; this does not add primal solves or refund failed work.
    pub fn evaluate_forward_with_constraints(&self, point: &[f64], control: &mut DesignControl,
        gate: &CancelGate) -> Result<ForwardDesignResponses, DesignError>
    {
        self.forward_query(point, control, gate, true)
    }

    fn forward_query(&self, point: &[f64], control: &mut DesignControl, gate: &CancelGate,
        include_constraints: bool) -> Result<ForwardDesignResponses, DesignError>
    {
        let prepared = self.prepare(point, control, gate)?;
        let mut cases = Vec::with_capacity(prepared.cases.len());
        let count = if include_constraints { self.constraints.len() } else { 0 };
        let mut rows = vec![None; count];
        for (index, case) in prepared.cases.iter().enumerate() {
            let solved = self.solve_case(&prepared, index, control, gate)?;
            let q: Vec<f64> = solved.network.models.iter().flat_map(|m| m.states())
                .map(|s| s.displacement_m_sqrt_kg).collect();
            let mut observations_m = Vec::with_capacity(case.targets.len());
            for target in &case.targets {
                checkpoint(gate)?;
                let map = &target.attachment;
                let start = solved.network.offsets[map.component];
                let component_q = &q[start..start + map.shapes.len()];
                // The same checked dot product and component-major order as
                // displacement_objective, with no objective residual arithmetic.
                observations_m.push(dot(&map.shapes, component_q).map_err(|e| case_error(index, e))?);
            }
            if include_constraints {
                for (row, constraint) in self.constraints.iter().enumerate() {
                    if constraint.case == index {
                        rows[row] = Some(responses::observe(constraint, &solved.network, &q,
                            &prepared.contacts, gate).map_err(|e| case_error(index, e))?);
                    }
                }
            }
            cases.push(ForwardDesignCase { observations_m, stored_energy_j: solved.stored_energy_j });
        }
        let constraints = rows.into_iter().map(|row|
            row.ok_or_else(|| bad("missing complete forward constraint row"))).collect::<Result<_,_>>()?;
        checkpoint(gate)?;
        Ok(ForwardDesignResponses {
            observations: ForwardDesignEvaluation { physical_parameters: prepared.parameters, cases,
                unassessed_response_constraints: self.constraints.len() - count },
            constraints,
        })
    }
}
