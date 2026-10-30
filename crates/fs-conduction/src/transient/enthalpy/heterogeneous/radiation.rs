//! Select an existing radiative boundary without changing the storage owner.

use super::{
    ConductionProblem, Cx, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyBackwardEuler, ThermalInterfaces,
};
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, PointwiseAmbientRadiation,
    PointwiseEnclosureRadiation,
    radiation::GrayDiffuseEnclosure,
    transient::enthalpy::{
        enclosure::EnthalpyEnclosureStepSolution, radiation::EnthalpyRadiationStepSolution,
    },
};

/// Explicitly chosen boundary approximation on the original enthalpy mesh.
///
/// The borrowed pointwise owners retain their prepared geometry/partition;
/// mean-patch variants retain the original caller's nonlinear exchange gates.
/// Every variant keeps the same conductivity, source, contact and caloric
/// inventory. No radiation law is inferred from a boundary or material name.
#[derive(Clone, Copy)]
pub enum EnthalpyRadiation<'a> {
    /// Only the original prescribed conduction boundary participates.
    None,
    /// Gray emission at each patch's area-mean temperature to a reservoir.
    Ambient {
        patches: &'a [AmbientRadiationPatch],
        config: AmbientRadiationConfig,
    },
    /// Closed reflected exchange at each patch's area-mean temperature.
    Enclosure {
        enclosure: &'a GrayDiffuseEnclosure,
        config: AmbientRadiationConfig,
    },
    /// Spatial P1 emission to the declared ambient reservoirs.
    PointwiseAmbient(&'a PointwiseAmbientRadiation<'a>),
    /// Spatial P1 emission and reflected, patch-uniform irradiation.
    PointwiseEnclosure(&'a PointwiseEnclosureRadiation<'a>),
}

/// Complete accepted result of the selected existing heat owner.
///
/// Mean-patch methods retain their applied inner boundary AND independent
/// physical endpoint audit. In those variants the inner Robin/Neumann heat
/// includes applied radiation; external consumers must use the separated
/// accounts in the original radiative report. Pointwise reports already live
/// in `EnthalpyStepSolution` and keep original convection separate.
#[derive(Debug, Clone)]
pub enum EnthalpyHeatStep {
    /// Prescribed-boundary or pointwise-radiative Newton endpoint.
    Native(EnthalpyStepSolution),
    /// Original complete mean-patch ambient endpoint and physical audit.
    Ambient(EnthalpyRadiationStepSolution),
    /// Original complete reflected mean-patch endpoint and physical audit.
    Enclosure(EnthalpyEnclosureStepSolution),
}

impl EnthalpyHeatStep {
    /// Original accepted enthalpy, phase, temperature and inner solve report.
    /// Consult the enclosing variant for a mean-patch physical energy audit.
    pub const fn conduction(&self) -> &EnthalpyStepSolution {
        match self {
            Self::Native(step) => step,
            Self::Ambient(step) => &step.conduction,
            Self::Enclosure(step) => &step.conduction,
        }
    }

    /// Complete physical storage-minus-net-input defect, in joules.
    pub const fn physical_energy_residual_j(&self) -> f64 {
        match self {
            Self::Native(step) => step.energy_residual_j,
            Self::Ambient(step) => step.physical_energy_residual_j,
            Self::Enclosure(step) => step.physical_energy_residual_j,
        }
    }
}

impl HeterogeneousEnthalpyBackwardEuler<'_, '_> {
    /// Evaluate the original physical residual at a supplied endpoint under
    /// the selected existing radiation owner. No history or nonlinear solve
    /// is advanced. An outer coupled solver must admit both the local balances
    /// and complete energy account before publishing its joint endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<super::super::EnthalpyEndpointEvaluation, EnthalpyError> {
        match radiation {
            EnthalpyRadiation::None => self.evaluate_endpoint(
                cx, problem, interfaces, old_specific_h, candidate_specific_h,
                dt_s, config, None,
            ),
            EnthalpyRadiation::Ambient { patches, config: radiation_config } => {
                radiation_config.validate()?;
                self.evaluate_endpoint(
                    cx, problem, interfaces, old_specific_h, candidate_specific_h,
                    dt_s, config, Some(patches),
                )
            }
            EnthalpyRadiation::Enclosure { enclosure, config: radiation_config } =>
                self.evaluate_endpoint_with_enclosure(
                    cx, problem, interfaces, old_specific_h, candidate_specific_h,
                    dt_s, config, enclosure, radiation_config,
                ),
            EnthalpyRadiation::PointwiseAmbient(owner) =>
                self.evaluate_endpoint_with_pointwise_ambient_radiation(
                    cx, problem, interfaces, old_specific_h, candidate_specific_h,
                    dt_s, config, owner,
                ),
            EnthalpyRadiation::PointwiseEnclosure(owner) =>
                self.evaluate_endpoint_with_pointwise_enclosure_radiation(
                    cx, problem, interfaces, old_specific_h, candidate_specific_h,
                    dt_s, config, owner,
                ),
        }
    }

    /// Advance the same immutable history through an explicitly selected
    /// existing radiation owner. All native residual, source-domain, energy,
    /// allocation, iteration and cancellation gates remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<EnthalpyHeatStep, EnthalpyError> {
        match radiation {
            EnthalpyRadiation::None => self
                .advance(cx, problem, interfaces, old_specific_h, dt_s, config)
                .map(EnthalpyHeatStep::Native),
            EnthalpyRadiation::Ambient {
                patches,
                config: radiation_config,
            } => self
                .advance_with_ambient_radiation(
                    cx,
                    problem,
                    interfaces,
                    old_specific_h,
                    dt_s,
                    config,
                    patches,
                    radiation_config,
                )
                .map(EnthalpyHeatStep::Ambient),
            EnthalpyRadiation::Enclosure {
                enclosure,
                config: radiation_config,
            } => self
                .advance_with_enclosure_radiation(
                    cx,
                    problem,
                    interfaces,
                    old_specific_h,
                    dt_s,
                    config,
                    enclosure,
                    radiation_config,
                )
                .map(EnthalpyHeatStep::Enclosure),
            EnthalpyRadiation::PointwiseAmbient(owner) => self
                .advance_with_pointwise_ambient_radiation(
                    cx,
                    problem,
                    interfaces,
                    old_specific_h,
                    dt_s,
                    config,
                    owner,
                )
                .map(EnthalpyHeatStep::Native),
            EnthalpyRadiation::PointwiseEnclosure(owner) => self
                .advance_with_pointwise_enclosure_radiation(
                    cx,
                    problem,
                    interfaces,
                    old_specific_h,
                    dt_s,
                    config,
                    owner,
                )
                .map(EnthalpyHeatStep::Native),
        }
    }
}
