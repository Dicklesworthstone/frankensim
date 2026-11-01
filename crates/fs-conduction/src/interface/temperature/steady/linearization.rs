//! Checked steady derivatives of the complete bulk and contact residual.

use fs_exec::Cx;
use crate::{ConductionError, ConductionProblem, EnergyBalance, LinearConfig};
use crate::adjoint::robin::{RobinGradient, RobinResponse};
use crate::interface::ContactResistanceGradient;
use crate::interface::temperature::tangent::ContactTemperatureLinearization;
use super::{TemperatureContacts, invalid, poll};

/// A revalidated steady endpoint and the complete nonsymmetric response.
/// No additional primal solve or derivative of solver iterations is performed.
pub struct ContactSteadyLinearization {
    response: RobinResponse,
    contact: ContactTemperatureLinearization,
    residual_norm_w: f64,
    energy: EnergyBalance,
    dirichlet_reactions_w: Vec<(usize, f64)>,
}

impl std::ops::Deref for ContactSteadyLinearization {
    type Target = RobinResponse;
    fn deref(&self) -> &Self::Target { &self.response }
}

/// One full transposed solve and its independent whole-law contact controls.
#[derive(Debug)]
pub struct ContactSteadyGradient {
    /// Source-load and selected Robin-reference/log(h) derivatives.
    pub transport: RobinGradient,
    /// dJ/ds for each named interface's R(T_A,T_B) -> s R(T_A,T_B), at s=1.
    pub resistance_multipliers: Vec<ContactResistanceGradient>,
}

impl ContactSteadyLinearization {
    /// Actual pinned contact states and both actions of their temperature tangent.
    #[must_use]
    pub const fn contact(&self) -> &ContactTemperatureLinearization { &self.contact }
    /// Euclidean free-node residual rechecked at the unchanged endpoint, W.
    #[must_use]
    pub const fn residual_norm_w(&self) -> f64 { self.residual_norm_w }
    /// Whole-domain source and reservoir balance at the same physical endpoint.
    #[must_use]
    pub const fn energy(&self) -> &EnergyBalance { &self.energy }
    /// Original Dirichlet reservoir heat entering at each fixed vertex, W.
    #[must_use]
    pub fn dirichlet_reactions_w(&self) -> &[(usize, f64)] { &self.dirichlet_reactions_w }

    /// Differentiate nodal/Robin observations through the full coupled residual.
    /// Contact controls use this SAME total adjoint, including both source
    /// temperature partials. Geometry, prescribed temperatures, source selection
    /// and every nonthermal context coordinate stay fixed. These controls have
    /// no additional physical time factor. No perturbed primal is executed.
    pub fn pullback_with_resistance_controls(
        &self, cx: &Cx<'_>, temperature_weights: &[f64], mean_weights: &[f64],
        heat_weights: &[f64], max_face_pairs: usize,
    ) -> Result<ContactSteadyGradient, ConductionError> {
        poll(cx, 0)?;
        if self.contact.trace_elements() > max_face_pairs {
            return Err(invalid("<steady-contact-adjoint>", "contact controls exceed the declared trace-element allowance"));
        }
        let transport = self.response.pullback(cx, temperature_weights, mean_weights, heat_weights)?;
        let resistance_multipliers = self.contact.physical_contacts().frozen_interfaces()
            .resistance_scale_pullback(cx, self.temperature(), &transport.nodal_load, max_face_pairs)?;
        Ok(ContactSteadyGradient { transport, resistance_multipliers })
    }
}

impl TemperatureContacts<'_, '_> {
    /// Revalidate and linearize an unchanged steady physical temperature field.
    ///
    /// The original bulk operator, source-backed contact law and prescribed
    /// values must satisfy the supplied absolute watt residual/heat targets.
    /// Those targets remain independent of the derivative solver's relative
    /// tolerance. An externally constructed result or a globally cancelling
    /// heat defect cannot bypass the physical endpoint checks.
    ///
    /// Matching laws retain face means; nonmatching laws retain pointwise
    /// common-refinement queries. The response includes bulk K'(T) and both
    /// contact temperature partials. Its tangent and true transpose use the
    /// existing bounded FGMRES response owner, even for a constant bulk tensor.
    /// Material/source corners do not acquire an invented unique derivative.
    /// Source uncertainty, shape changes, pressure evolution and prescribed-
    /// temperature controls are outside this local constitutive derivative.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_steady_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
        linear: LinearConfig, residual_tolerance_w: f64, energy_tolerance_w: f64,
        regions: &[&str], max_face_pairs: usize,
    ) -> Result<ContactSteadyLinearization, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        if !(residual_tolerance_w.is_finite() && residual_tolerance_w > 0.0)
            || !(energy_tolerance_w.is_finite() && energy_tolerance_w > 0.0)
        {
            return Err(invalid("<steady-contact-adjoint>", "finite positive watt residual and heat targets are required"));
        }
        let contact = self.linearize_temperature(cx, problem.boundary, temperature, max_face_pairs)?;
        let (response, residual_norm_w, energy, dirichlet_reactions_w) =
            RobinResponse::for_steady_contact(cx, problem, &contact, linear,
                residual_tolerance_w, energy_tolerance_w, regions)?;
        Ok(ContactSteadyLinearization { response, contact, residual_norm_w,
            energy, dirichlet_reactions_w })
    }
}
