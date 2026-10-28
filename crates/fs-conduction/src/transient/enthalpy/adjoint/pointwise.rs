//! Explicit adjoints of accepted spatial P1 radiation endpoints.
//!
//! The caller supplies the original nonlinear residual target. Binding
//! reevaluates the complete physical load and energy balance, then retains
//! exactly the local emission blocks and reflected patch response used by
//! the pointwise Newton owner. No shifted Robin/flux surrogate is bound.
//!
//! Emissivity remains a frozen scalar declaration on each trace. State
//! derivatives include T_h^4 quadrature and enclosure reflection, not an
//! undeclared temperature derivative of the emissivity card query.

use super::{
    EnthalpyAdjointError, EnthalpyBackwardEuler, EnthalpyStepConfig,
    EnthalpyStepLinearization, EnthalpyStepSolution, PointwiseBoundary,
};
use crate::{
    ConductionProblem, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    ThermalInterfaces,
};
use fs_exec::Cx;

/// Conservative additional scalar/index-word allowance for a retained
/// pointwise endpoint and its transpose/parameter work vectors.
///
/// Counts include the local face blocks, sparse patch geometry, small dense
/// reflected response and temporary nodal/patch vectors. They exclude source
/// strings and are not an operating-system RSS guarantee. Checked arithmetic
/// precedes tangent allocation; callers convert entries to bytes separately.
pub fn workspace_entries(
    vertices: usize,
    boundary_faces: usize,
    patches: usize,
) -> Result<usize, EnthalpyAdjointError> {
    if vertices == 0 || boundary_faces == 0 || !(1..=64).contains(&patches) {
        return Err(EnthalpyAdjointError::InvalidInput(
            "pointwise adjoint requires nonempty geometry and 1..=64 radiation patches",
        ));
    }
    let overflow = || EnthalpyAdjointError::InvalidInput("pointwise adjoint workspace overflow");
    let nodal = vertices.checked_mul(8).ok_or_else(overflow)?;
    let faces = boundary_faces.checked_mul(48).ok_or_else(overflow)?;
    let matrix = patches.checked_mul(patches).and_then(|n| n.checked_mul(4))
        .ok_or_else(overflow)?;
    let patch = patches.checked_mul(32).ok_or_else(overflow)?;
    nodal.checked_add(faces).and_then(|n| n.checked_add(matrix))
        .and_then(|n| n.checked_add(patch)).ok_or_else(overflow)
}

fn admit(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    patches: usize,
    physical_target_j: f64,
    max_radiation_entries: usize,
) -> Result<(), EnthalpyAdjointError> {
    super::poll(cx, 0)?;
    if !physical_target_j.is_finite() || physical_target_j < 0.0 {
        return Err(EnthalpyAdjointError::InvalidInput(
            "finite nonnegative original pointwise residual target required",
        ));
    }
    let entries = workspace_entries(
        problem.mesh.vertex_count(), problem.mesh.boundary().len(), patches,
    )?;
    if entries > max_radiation_entries {
        return Err(EnthalpyAdjointError::InvalidInput(
            "pointwise radiation adjoint workspace exceeds the declared entry cap",
        ));
    }
    Ok(())
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Bind a physically accepted pointwise ambient-radiation endpoint.
    ///
    /// The owner must borrow this exact mesh and current physical boundary.
    /// Every radiating node must be strictly inside retained emissivity
    /// temperature support; chart/material derivatives must likewise be
    /// single-valued. The unchanged original Newton target may be zero.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
        max_radiation_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        admit(cx, problem, radiation.patch_count(), physical_target_j, max_radiation_entries)?;
        if !radiation.matches(problem.mesh, problem.boundary)
            || accepted.pointwise_radiation.is_none()
            || accepted.pointwise_enclosure_radiation.is_some()
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "ambient pointwise owner, boundary and accepted radiation model must agree",
            ));
        }
        self.linearize_accepted_model_with_target(
            cx, problem, interfaces, old_h, dt_s, config, accepted,
            Some(physical_target_j), Some(PointwiseBoundary::Ambient(radiation)),
        )
    }

    /// Bind the complete reflected pointwise enthalpy Jacobian at an accepted
    /// endpoint, keeping original Robin/Neumann laws and all bath columns.
    /// The returned base transport owner supports the existing material,
    /// contact, convective-feedback, storage and source contractions.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseEnclosureRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
        max_radiation_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        admit(cx, problem, radiation.patch_count(), physical_target_j, max_radiation_entries)?;
        if !radiation.matches(problem.mesh, problem.boundary)
            || accepted.pointwise_enclosure_radiation.is_none()
            || accepted.pointwise_radiation.is_some()
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "enclosure pointwise owner, boundary and accepted radiation model must agree",
            ));
        }
        self.linearize_accepted_model_with_target(
            cx, problem, interfaces, old_h, dt_s, config, accepted,
            Some(physical_target_j), Some(PointwiseBoundary::Enclosure(radiation)),
        )
    }
}
