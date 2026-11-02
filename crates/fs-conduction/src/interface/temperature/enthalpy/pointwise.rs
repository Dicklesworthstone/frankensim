//! Accepted contact/radiation derivatives use the same physical enthalpy tangent.
//!
//! These bindings compose the original pointwise radiation and source-backed
//! contact owners. They do not differentiate Newton iterations or advance the
//! primal history. Radiation and contact contributions occur exactly once.

use super::{ContactEnthalpyGradient, ContactEnthalpyLinearization, ContactEnthalpyStep};
use crate::interface::temperature::{
    TemperatureContacts, tangent::ContactTemperatureLinearization,
};
use crate::transient::enthalpy::{
    EnthalpyStepConfig, adjoint::EnthalpyAdjointError,
    heterogeneous::HeterogeneousEnthalpyBackwardEuler,
};
use crate::{
    ConductionProblem, LinearConfig, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
};
use fs_exec::Cx;

/// Physical thermal, whole-contact-law and pointwise radiation controls from
/// one transpose solve. Radiation controls hold geometry and view factors fixed.
#[derive(Debug)]
pub struct ContactEnthalpyPointwiseGradient {
    /// History, heat sources, reference density and resistance multipliers.
    pub thermal: ContactEnthalpyGradient,
    /// Original emissivity and reservoir-temperature derivatives in retained
    /// trace order. No inferred emissivity state law is added.
    pub radiation: crate::radiation::pointwise::PointwiseRadiationParameterGradient,
}

impl TemperatureContacts<'_, '_> {
    pub(super) fn rebind_pointwise_enthalpy_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        accepted: &ContactEnthalpyStep,
        physical_target_j: f64,
        max_face_pairs: usize,
    ) -> Result<ContactTemperatureLinearization, EnthalpyAdjointError> {
        super::poll(cx, 0)?;
        if accepted.radiation.is_some()
            || accepted.combined_boundary.is_some()
            || !physical_target_j.is_finite()
            || physical_target_j <= 0.0
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "pointwise contact adjoint requires its original positive physical target and no mean-patch boundary",
            ));
        }
        let contact = self.linearize_temperature(
            cx,
            problem.boundary,
            &accepted.conduction.temperature,
            max_face_pairs,
        )?;
        if contact.physical_contacts().states() != accepted.physical_contacts.states()
            || contact.physical_contacts().quadrature_states()
                != accepted.physical_contacts.quadrature_states()
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "rebound pointwise contact source state differs from the accepted physical endpoint",
            ));
        }
        Ok(contact)
    }

    /// Recheck an accepted local-P1 ambient-radiation/contact endpoint and bind
    /// its complete Jacobian, including dR_contact/dT and dT/dh on latent states.
    ///
    /// The original physical residual target remains authoritative. The
    /// radiation owner must match this exact mesh and boundary. Geometry,
    /// contact pressure and declared emissivities are frozen; the pointwise
    /// owner retains all admitted emission derivatives. The existing chart
    /// derivative owner refuses nonsmooth enthalpy knots.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy_with_pointwise_ambient_radiation<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_target_j: f64,
        radiation: &PointwiseAmbientRadiation<'_>,
        accepted: ContactEnthalpyStep,
        max_radiation_entries: usize,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyLinearization<'m>, EnthalpyAdjointError> {
        let contact = self.rebind_pointwise_enthalpy_contacts(
            cx,
            problem,
            &accepted,
            physical_target_j,
            max_face_pairs,
        )?;
        let step = engine
            .linearize_accepted_with_pointwise_ambient_radiation(
                cx,
                problem,
                Some(contact.physical_contacts().frozen_interfaces()),
                old_h,
                dt_s,
                config,
                radiation,
                accepted.conduction,
                physical_target_j,
                max_radiation_entries,
            )?
            .with_contact_temperature_feedback(cx, &contact)?;
        Ok(ContactEnthalpyLinearization { step, contact })
    }

    /// Bind the same accepted enthalpy/contact model with full pointwise
    /// enclosure reflection, retaining its original reservoir and spectral
    /// columns where the supplied enclosure owner admits them.
    ///
    /// The same endpoint, smoothness and work requirements as the ambient
    /// binding apply. The reflected response stays matrix-free; this does not
    /// expand a dense nodal radiation matrix or run another primal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy_with_pointwise_enclosure_radiation<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_target_j: f64,
        radiation: &PointwiseEnclosureRadiation<'_>,
        accepted: ContactEnthalpyStep,
        max_radiation_entries: usize,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyLinearization<'m>, EnthalpyAdjointError> {
        let contact = self.rebind_pointwise_enthalpy_contacts(
            cx,
            problem,
            &accepted,
            physical_target_j,
            max_face_pairs,
        )?;
        let step = engine
            .linearize_accepted_with_pointwise_enclosure_radiation(
                cx,
                problem,
                Some(contact.physical_contacts().frozen_interfaces()),
                old_h,
                dt_s,
                config,
                radiation,
                accepted.conduction,
                physical_target_j,
                max_radiation_entries,
            )?
            .with_contact_temperature_feedback(cx, &contact)?;
        Ok(ContactEnthalpyLinearization { step, contact })
    }
}

impl ContactEnthalpyLinearization<'_> {
    /// Differentiate one objective with respect to thermal history, sources,
    /// whole resistance laws and radiation controls using one complete transpose.
    /// Convert temperature objective weights with `temperature_pullback` first.
    ///
    /// Radiation leaves the solid, so its control contraction is the negative
    /// of the incoming nodal-load contraction. Internal contact heat has no
    /// additional external source account. Direct objective dependence on a
    /// control remains the caller's responsibility.
    pub fn pullback_with_pointwise_radiation_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: LinearConfig,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyPointwiseGradient, EnthalpyAdjointError> {
        let transport = self.step.transport();
        if transport.primal().pointwise_radiation.is_none()
            && transport.primal().pointwise_enclosure_radiation.is_none()
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "pointwise radiation controls require a retained pointwise contact endpoint",
            ));
        }
        let thermal = self.pullback_with_resistance_controls(cx, seed, linear, max_face_pairs)?;
        let mut outward_weights =
            transport.nodal_load_pullback(cx, &thermal.storage.transport.adjoint)?;
        for (vertex, weight) in outward_weights.iter_mut().enumerate() {
            if vertex % crate::assemble::ASSEMBLY_TILE == 0 {
                super::poll(cx, vertex)?;
            }
            *weight = -*weight;
        }
        let radiation = transport.pointwise_parameter_pullback(cx, &outward_weights)?;
        Ok(ContactEnthalpyPointwiseGradient { thermal, radiation })
    }
}
