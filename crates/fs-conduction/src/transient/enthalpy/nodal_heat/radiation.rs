//! Simultaneous native pointwise radiation and held signed nodal interval heat.
//!
//! This selects existing radiation owners without changing their geometry,
//! caloric identities, source claims, quadrature, reflection or exact tangent.
//! The nodal source is constant during Newton; a coupled owner recomputes its
//! physical source and uses endpoint evaluation before admitting a trajectory.
//! Pointwise radiation and convection keep separate accounts. Mean-patch
//! frozen-boundary drivers are different discrete models and are not selected
//! by this API. The original backward-Euler order and fixed geometry remain.

use super::*;
use super::super::heterogeneous::EnthalpyRadiation;

impl EnthalpyBackwardEuler<'_, '_> {
    /// Solve signed nodal heat together with original pointwise ambient or
    /// reflected radiation, including material/phase-dependent emissivity.
    ///
    /// The selected owner must retain this exact mesh and boundary; material
    /// owners must also retain every original nodal caloric object. Local
    /// source residuals, prescribed-node reactions and the independent global
    /// balance all include radiation at the actual returned enthalpy.
    /// The original native Newton tangent supplies emission and reflection
    /// derivatives; no radiation fixed point or duplicate boundary is added.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat_and_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<NodalHeatStep, EnthalpyError> {
        let radiation = pointwise(self, cx, problem, radiation)?;
        self.advance_nodal_heat(
            cx, problem, interfaces, old_specific_h, dt_s, config, nodal_heat_j, radiation,
        )
    }

    /// Re-evaluate a supplied endpoint with its freshly recomputed signed
    /// source and original pointwise radiation. This performs no time step.
    ///
    /// The calling physical owner must admit both local and global residuals;
    /// equal-and-opposite source errors cannot pass on a total balance alone.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat_and_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<NodalHeatEndpoint, EnthalpyError> {
        let radiation = pointwise(self, cx, problem, radiation)?;
        self.evaluate_nodal_heat(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, nodal_heat_j, radiation,
        )
    }
}

pub(super) fn pointwise<'a>(
    storage: &EnthalpyBackwardEuler<'_, '_>,
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    radiation: EnthalpyRadiation<'a>,
) -> Result<Option<PointwiseBoundary<'a>>, EnthalpyError> {
    poll(cx, 0)?;
    let selected = match radiation {
        EnthalpyRadiation::None => return Ok(None),
        EnthalpyRadiation::PointwiseAmbient(owner) => {
            if !owner.matches(problem.mesh, problem.boundary) {
                return Err(EnthalpyError::InvalidInput(
                    "nodal heat radiation requires the original mesh and boundary",
                ));
            }
            PointwiseBoundary::Ambient(owner)
        }
        EnthalpyRadiation::PointwiseEnclosure(owner) => {
            if !owner.matches(problem.mesh, problem.boundary) {
                return Err(EnthalpyError::InvalidInput(
                    "nodal heat enclosure requires the original mesh and boundary",
                ));
            }
            PointwiseBoundary::Enclosure(owner)
        }
        EnthalpyRadiation::MaterialAmbient(owner) => {
            if !std::ptr::eq(storage.mesh, problem.mesh)
                || !owner.matches(problem.mesh, problem.boundary)
                || !owner.matches_curves(cx, |v| storage.curve_for_vertex(v))?
            {
                return Err(EnthalpyError::InvalidInput(
                    "nodal heat material radiation requires the original mesh, boundary and nodal caloric objects",
                ));
            }
            PointwiseBoundary::MaterialAmbient(owner)
        }
        EnthalpyRadiation::MaterialEnclosure(owner) => {
            if !std::ptr::eq(storage.mesh, problem.mesh)
                || !owner.matches(problem.mesh, problem.boundary)
                || !owner.matches_curves(cx, |v| storage.curve_for_vertex(v))?
            {
                return Err(EnthalpyError::InvalidInput(
                    "nodal heat material enclosure requires the original mesh, boundary and nodal caloric objects",
                ));
            }
            PointwiseBoundary::MaterialEnclosure(owner)
        }
        EnthalpyRadiation::Ambient { .. } | EnthalpyRadiation::Enclosure { .. } => {
            return Err(EnthalpyError::InvalidInput(
                "signed nodal heat requires an original pointwise radiation owner; mean-patch radiation is a different discrete model",
            ));
        }
    };
    Ok(Some(selected))
}
