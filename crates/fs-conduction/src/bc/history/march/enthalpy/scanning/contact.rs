//! Retained material-contact feedback for the original moving-source trajectory.
//!
//! The scan owns source projection and publication. TemperatureContacts owns all
//! interface geometry, source queries, nonlinear coupling and physical audits;
//! no frozen contact snapshot is substituted for that owner.

use super::{ScanningEnthalpyRun, conduction, invalid, poll};
use crate::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionProblem};
use crate::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use crate::interface::temperature::enthalpy::{ContactEnthalpyNewtonConfig, ContactEnthalpyStep};
use crate::transient::enthalpy::{EnthalpyError, EnthalpyStepConfig};
use crate::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError,
};
use fs_exec::Cx;

#[derive(Debug, Clone, Copy)]
pub(super) enum ScanContactCoupling {
    Legacy(ContactIterationConfig),
    Newton(ContactEnthalpyNewtonConfig),
}

pub(super) struct ScanContactBinding<'m, 'c> {
    pub(super) contacts: TemperatureContacts<'m, 'c>,
    pub(super) coupling: ScanContactCoupling,
    pub(super) max_face_pairs: usize,
}

impl ScanContactBinding<'_, '_> {
    pub(super) fn physical_residual_tolerance_j(&self) -> f64 {
        match self.coupling {
            ScanContactCoupling::Legacy(config) => config.residual_tolerance_j,
            ScanContactCoupling::Newton(config) => config.residual_tolerance_j,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance(
        &self, cx: &Cx<'_>, storage: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>, old: &[f64], dt_s: f64, step: EnthalpyStepConfig,
        radiation: Option<(&[AmbientRadiationPatch], AmbientRadiationConfig)>,
    ) -> Result<ContactEnthalpyStep, HeterogeneousEnthalpyError> {
        match self.coupling {
            ScanContactCoupling::Newton(config) => {
                if config.step != step {
                    return Err(invalid("scan advance controls must match the retained contact enthalpy Newton policy"));
                }
                if radiation.is_some() {
                    return Err(invalid("coupled Newton contacts require pointwise radiation, not the scan's area-mean radiation"));
                }
                Ok(self.contacts.advance_enthalpy_newton(cx, storage, problem, old, dt_s, config)?)
            }
            ScanContactCoupling::Legacy(coupling) => match radiation {
                None => Ok(self.contacts.advance_enthalpy(
                    cx, storage, problem, old, dt_s, step, coupling)?),
                Some((patches, config)) => Ok(self.contacts.advance_enthalpy_with_ambient_radiation(
                    cx, storage, problem, old, dt_s, step, coupling, patches, config)?),
            },
        }
    }
}

impl<'h, 'm, 'c> ScanningEnthalpyRun<'h, 'm, 'c> {
    /// Retain the complete fixed and source-backed contact set for this scan.
    /// Moving deposition, background heating, original phase charts, boundary
    /// histories and optional ambient radiation participate in each contact trial.
    /// The source owner updates its temperature coordinates; its geometry,
    /// pressure and other explicitly frozen source coordinates stay fixed.
    ///
    /// Attach once before the first accepted step. Put every fixed interface in
    /// the supplied contact set; a separately declared base set is refused.
    /// The positive face-pair cap bounds later contact tangent preparation.
    /// Geometry admission retains the original contact-construction allowance.
    /// Material-phase contact feedback requires with_contacts_newton instead.
    ///
    /// Construction checks the initial source queries without advancing time.
    pub fn with_contacts(
        mut self, cx: &Cx<'_>, contacts: TemperatureContacts<'m, 'c>,
        coupling: ContactIterationConfig, max_face_pairs: usize,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        self.admit_contacts(cx, &contacts, max_face_pairs)?;
        if coupling.max_iterations == 0 || !coupling.relaxation.is_finite()
            || coupling.relaxation <= 0.0 || coupling.relaxation > 1.0
            || !coupling.residual_tolerance_j.is_finite() || coupling.residual_tolerance_j <= 0.0
        {
            return Err(invalid("scanning contact needs positive iteration/residual budgets and relaxation in (0,1]"));
        }
        let boundary = self.run.problem.boundary.at(cx, self.run.time_s()).map_err(conduction)?;
        contacts.evaluate(cx, &boundary, &self.run.temperature).map_err(conduction)?;
        poll(cx, 0).map_err(conduction)?;
        self.contacts = Some(ScanContactBinding {
            contacts, coupling: ScanContactCoupling::Legacy(coupling), max_face_pairs,
        });
        Ok(self)
    }

    /// Couple moving absorbed power to the original full enthalpy/contact Newton
    /// residual, including independent material liquid-fraction coordinates.
    /// Both bodies keep their own caloric charts; latent evolution remains active
    /// when a temperature derivative vanishes. Every accepted report retains the
    /// source-backed contact states, counted work and independent joule audit.
    ///
    /// The complete policy and its face-pair cap remain fixed for continuation.
    /// Advance calls must pass exactly config.step. Area-mean ambient radiation
    /// is refused because the coupled Newton owner requires pointwise radiation.
    /// No thermal state, physical time or step allowance changes on attachment.
    pub fn with_contacts_newton(
        mut self, cx: &Cx<'_>, contacts: TemperatureContacts<'m, 'c>,
        config: ContactEnthalpyNewtonConfig,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        self.admit_contacts(cx, &contacts, config.max_face_pairs)?;
        if self.radiation.is_some() {
            return Err(invalid("coupled Newton scan contacts require pointwise radiation; area-mean scan radiation is incompatible"));
        }
        config.validate()?;
        config.step.newton.validate().map_err(EnthalpyError::Newton)?;
        if config.step.max_newton_iterations == 0 || !config.step.energy_tolerance_j.is_finite()
            || config.step.energy_tolerance_j <= 0.0
        {
            return Err(invalid("scanning contact requires positive nonlinear work and finite positive energy tolerance"));
        }
        let boundary = self.run.problem.boundary.at(cx, self.run.time_s()).map_err(conduction)?;
        contacts.evaluate_enthalpy(cx, &boundary, &self.run.storage, &self.run.specific_h)?;
        poll(cx, 0).map_err(conduction)?;
        self.contacts = Some(ScanContactBinding {
            contacts, coupling: ScanContactCoupling::Newton(config),
            max_face_pairs: config.max_face_pairs,
        });
        Ok(self)
    }

    /// Original complete contact model, never a frozen endpoint replacement.
    #[must_use]
    pub fn contacts(&self) -> Option<&TemperatureContacts<'m, 'c>> {
        self.contacts.as_ref().map(|binding| &binding.contacts)
    }

    /// Retained temperature-only iteration policy; absent for coupled Newton.
    #[must_use]
    pub fn contact_coupling(&self) -> Option<ContactIterationConfig> {
        match self.contacts.as_ref().map(|binding| binding.coupling) {
            Some(ScanContactCoupling::Legacy(config)) => Some(config),
            _ => None,
        }
    }

    /// Complete original coupled Newton policy, when selected.
    #[must_use]
    pub fn contact_newton_config(&self) -> Option<ContactEnthalpyNewtonConfig> {
        match self.contacts.as_ref().map(|binding| binding.coupling) {
            Some(ScanContactCoupling::Newton(config)) => Some(config),
            _ => None,
        }
    }

    fn admit_contacts(&self, cx: &Cx<'_>, contacts: &TemperatureContacts<'m, 'c>,
        max_face_pairs: usize) -> Result<(), HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() != 0 || self.contacts.is_some()
            || self.run.problem.interfaces.is_some() || max_face_pairs == 0
        {
            return Err(invalid("attach the complete contact set once to an unadvanced scan with no separate interfaces and a positive tangent face cap"));
        }
        contacts.check_mesh(self.run.problem.boundary.mesh()).map_err(conduction)?;
        Ok(())
    }
}
