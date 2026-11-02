//! Original material emissivity controls on the complete contact/h transpose.
use super::*;
use crate::LinearConfig;
use crate::transient::enthalpy::adjoint::EnthalpyAdjointError;

#[derive(Debug)]
pub struct ContactEnthalpyMaterialAmbientGradient {
    pub thermal: ContactEnthalpyGradient,
    pub radiation:
        crate::radiation::pointwise::material::adjoint::MaterialRadiationParameterGradient,
}

#[derive(Debug)]
pub struct ContactEnthalpyMaterialEnclosureGradient {
    pub thermal: ContactEnthalpyGradient,
    pub radiation:
        crate::radiation::pointwise_enclosure::material::adjoint::MaterialEnclosureParameterGradient,
}

impl TemperatureContacts<'_, '_> {
    /// Rebind material enclosure reflection and the complete original contact
    /// law at one accepted endpoint; both owners require smooth source branches.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy_with_material_enclosure_radiation<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_target_j: f64,
        radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
        accepted: ContactEnthalpyStep,
        max_radiation_entries: usize,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyLinearization<'m>, EnthalpyAdjointError> {
        let contact = self.rebind_pointwise_enthalpy_contacts(cx, problem, engine,
            &accepted, physical_target_j, max_face_pairs)?;
        let step = engine.linearize_accepted_with_material_enclosure_radiation(
            cx, problem, Some(contact.physical_contacts().frozen_interfaces()), old_h, dt_s,
            config, radiation, accepted.conduction, physical_target_j, max_radiation_entries)?
            .with_contact_temperature_feedback(cx, &contact)?;
        Ok(ContactEnthalpyLinearization { step, contact })
    }

    /// Recheck and bind source-dependent ambient emissivity together with the
    /// complete contact law. Both original source owners enforce strict
    /// two-sided admission; the forward generalized slope is not an adjoint.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy_with_material_ambient_radiation<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_target_j: f64,
        radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
        accepted: ContactEnthalpyStep,
        max_radiation_entries: usize,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyLinearization<'m>, EnthalpyAdjointError> {
        let contact = self.rebind_pointwise_enthalpy_contacts(
            cx,
            problem,
            engine,
            &accepted,
            physical_target_j,
            max_face_pairs,
        )?;
        let step = engine
            .linearize_accepted_with_material_ambient_radiation(
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
    /// Return storage, contact resistance, original enclosure emissivity and
    /// open-bath controls from one physical transpose solve.
    pub fn pullback_with_material_enclosure_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: LinearConfig,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyMaterialEnclosureGradient, EnthalpyAdjointError> {
        let thermal = self.pullback_with_resistance_controls(cx, seed, linear, max_face_pairs)?;
        let transport = self.step.transport();
        let mut weights = transport.nodal_load_pullback(cx, &thermal.storage.transport.adjoint)?;
        for (vertex, weight) in weights.iter_mut().enumerate() {
            if vertex % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            *weight = -*weight;
        }
        let radiation = transport.material_enclosure_parameter_pullback(cx, &weights)?;
        Ok(ContactEnthalpyMaterialEnclosureGradient { thermal, radiation })
    }

    /// Return thermal, density, contact-resistance, material-emissivity and bath
    /// controls from one physical transpose solve. Direct objective control
    /// dependence remains with the caller.
    pub fn pullback_with_material_ambient_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: LinearConfig,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyMaterialAmbientGradient, EnthalpyAdjointError> {
        let thermal = self.pullback_with_resistance_controls(cx, seed, linear, max_face_pairs)?;
        let transport = self.step.transport();
        let mut weights = transport.nodal_load_pullback(cx, &thermal.storage.transport.adjoint)?;
        for (vertex, weight) in weights.iter_mut().enumerate() {
            if vertex % crate::assemble::ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            *weight = -*weight;
        }
        let radiation = transport.material_ambient_parameter_pullback(cx, &weights)?;
        Ok(ContactEnthalpyMaterialAmbientGradient { thermal, radiation })
    }
}
