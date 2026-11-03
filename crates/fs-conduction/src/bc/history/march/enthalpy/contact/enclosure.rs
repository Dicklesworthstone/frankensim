//! Original material contact and reflected radiation share one h history.

use super::*;
use crate::{
    ConductionProblem, MaterialEnclosureRadiationPatch, MaterialPointwiseEnclosureRadiation,
    PointwiseEnclosureRadiation, SurfaceEmissivity,
};
use super::super::radiation::enclosure::MaterialEnclosureHistoryBinding;
use fs_scenario::signal::TimeSignal;

impl<'h, 'm, 'c, 'i> ContactEnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    /// Attach original phase-dependent emissivity and closed/open reflection to
    /// the coupled enthalpy Newton history. Each candidate solves contact R,
    /// source k, material emission, reflection and reference-mass storage together.
    ///
    /// An open enclosure requires one independent black-bath temperature
    /// schedule; a closed enclosure requires `None`. Initial storage/source
    /// identity and partition admission are shared with the enclosure-only
    /// history. No initial contact or radiation query advances physical time.
    ///
    /// # Errors
    /// Non-Newton coupling, advanced history, duplicate radiative ownership,
    /// original mesh/chart/source/bath refusal or cancellation.
    pub fn with_material_enclosure(
        mut self, cx: &Cx<'_>, base: &'h PointwiseEnclosureRadiation<'h>,
        patches: Vec<MaterialEnclosureRadiationPatch<'c>>,
        ambient_temperature: Option<TimeSignal>,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if self.newton_config().is_none() || self.radiation.is_some()
            || self.material_enclosure.is_some()
        {
            return Err(invalid("material enclosure requires coupled enthalpy Newton and exclusive radiation ownership"));
        }
        self.material_enclosure = Some(MaterialEnclosureHistoryBinding::new(
            cx, &self.run, base, patches, ambient_temperature,
        )?);
        Ok(self)
    }

    /// Original pinned sources in the original enclosure's surface order.
    /// `None` means that no material enclosure was attached.
    #[must_use]
    pub fn material_enclosure_patches(&self) -> Option<&[MaterialEnclosureRadiationPatch<'c>]> {
        self.material_enclosure.as_ref().map(MaterialEnclosureHistoryBinding::patches)
    }

    /// Preparation-time emissivity answers; actual accepted nodal source states
    /// remain in each step's material pointwise enclosure report.
    #[must_use]
    pub fn material_enclosure_reference_emissivities(&self) -> Option<&[SurfaceEmissivity]> {
        self.material_enclosure.as_ref().map(MaterialEnclosureHistoryBinding::reference_emissivities)
    }

    /// Sole black-bath history for an open material enclosure. A closed enclosure
    /// has no external radiative reservoir and returns `None`.
    #[must_use]
    pub fn material_enclosure_ambient_temperature(&self) -> Option<&TimeSignal> {
        self.material_enclosure.as_ref().and_then(MaterialEnclosureHistoryBinding::ambient_temperature)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve_material_enclosure(
        &self, cx: &Cx<'_>, radiation: &MaterialEnclosureHistoryBinding<'h, 'c>,
        problem: ConductionProblem<'_>, old: &[f64], dt_s: f64, time_s: f64,
        step: EnthalpyStepConfig, before: bool,
    ) -> Result<ContactEnthalpyStep, HeterogeneousEnthalpyError> {
        let config = self.newton_config_for(step)?;
        let base = radiation.sample(cx, problem.mesh, problem.boundary, time_s, before)?;
        let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, radiation.patches())
            .map_err(conduction)?;
        Ok(self.contacts.advance_enthalpy_newton_with_material_enclosure_radiation(
            cx, &self.run.storage, problem, old, dt_s, config, &owner,
        )?)
    }
}
