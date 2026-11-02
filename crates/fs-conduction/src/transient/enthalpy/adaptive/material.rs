//! Adaptive material emissivity through the original enthalpy time controller.

use std::ops::ControlFlow;

use fs_exec::Cx;

use super::{
    AdaptiveEnthalpyReport, AdaptiveEnthalpyState, EnthalpyAcceptedTrial, EnthalpyError,
};
use crate::{ConductionProblem, MaterialPointwiseAmbientRadiation,
    MaterialPointwiseEnclosureRadiation, ThermalInterfaces};

impl AdaptiveEnthalpyState<'_, '_, '_> {
    /// Adapt material-dependent local emission and reflected enclosure exchange
    /// with the original h, T and phase error controller. Every coarse and fine
    /// endpoint resolves emissivity and refactors the actual reflection system.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_material_enclosure_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &MaterialPointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_material_enclosure_radiation_observed(
            cx, problem, interfaces, radiation, end_time_s, max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe only complete accepted fine pairs. Internal surface exchange
    /// cancels in the physical account; an explicit bath contributes its signed
    /// heat once. Rejected comparisons retain neither heat nor material history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_material_enclosure_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &MaterialPointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        self.advance_with(
            cx, end_time_s, max_trials,
            |h, dt| storage.advance_with_material_enclosure_radiation(
                cx, problem, interfaces, h, dt, config, radiation,
            ),
            observe,
        )
    }

    /// Adapt local P1 radiation with emissivity resolved from each node's
    /// original temperature and liquid fraction. The shared controller checks
    /// h, T and phase error independently, including isothermal latent storage.
    /// Source claims, geometry and reservoirs stay fixed for this event window;
    /// the material state is refreshed at every nonlinear endpoint candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_material_ambient_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &MaterialPointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_material_ambient_radiation_observed(
            cx,
            problem,
            interfaces,
            radiation,
            end_time_s,
            max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe complete accepted fine pairs with their actual nodal emissivity
    /// receipts. Comparison and rejected heat never enter accepted history.
    /// Stop/resume and event clipping retain the original controller semantics.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_material_ambient_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &MaterialPointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        self.advance_with(
            cx,
            end_time_s,
            max_trials,
            |h, dt| {
                storage.advance_with_material_ambient_radiation(
                    cx, problem, interfaces, h, dt, config, radiation,
                )
            },
            observe,
        )
    }
}
