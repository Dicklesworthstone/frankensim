//! Accepted capacity endpoints with the full source-backed R(T) Jacobian.
use super::*;
use crate::interface::temperature::{TemperatureContacts, tangent::ContactCapacityLinearization};

impl BackwardEuler<'_> {
    /// Revalidate and differentiate an unchanged R(T)-coupled capacity endpoint.
    ///
    /// Uses the original source pins and the caller's fixed absolute residual
    /// and energy gates. The candidate is never replaced or re-solved. The
    /// complete nonsymmetric Jacobian and its transpose use the existing
    /// bounded FGMRES response owner, including nonlinear bulk conductivity.
    ///
    /// Contact pressure, geometry, source selection and prescribed exterior
    /// laws remain fixed. A kink or one-sided source domain refuses even for
    /// a zero goal. `nodal_load` is already in watts: history/source/capacity
    /// and resistance-scale contractions must not multiply it by dt again.
    /// The face cap includes fixed contacts; it is additional to mesh/solver
    /// admission, not a total peak-memory limit.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_temperature_contact_endpoint<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>, old: &[f64], candidate: &[f64],
        dt_s: f64, config: StepConfig, residual_tolerance_j: f64,
        regions: &[&str], max_face_pairs: usize,
    ) -> Result<ContactCapacityLinearization<'a>, ConductionError> {
        poll(cx, 0)?;
        contacts.check_mesh(problem.mesh)?;
        let contact = contacts.linearize_temperature(cx, problem.boundary, candidate, max_face_pairs)?;
        let mut step = self.linearize_prescribed_endpoint(
            cx, problem, Some(contact.physical_contacts().frozen_interfaces()), old,
            candidate, dt_s, config, residual_tolerance_j, regions,
        )?;
        step.response = step.response.with_contact_temperature_feedback(cx, &contact)?;
        poll(cx, 0)?;
        Ok(ContactCapacityLinearization::new(step, contact))
    }
}
