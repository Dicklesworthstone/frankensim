//! Thermal numerical reactions through the original embedded-support owner.
//!
//! The original Nitsche functional uses C eps(u) n - gamma (u-g). Prescribed
//! free strain subtracts C eps_free n on the SAME selected support quadrature.
//! Its value and direct material-scale derivative change; its displacement
//! derivative remains the original lifting transpose.

use std::ops::ControlFlow;

use super::{PreparedEigenstrain3, poll};
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::dirichlet::EmbeddedReaction3;
use crate::elastic3::{CutElasticity3, ElasticityError3};

impl CutElasticity3 {
    /// Nitsche-consistent support reaction from the original thermal state.
    ///
    /// Integrates `mode dot [scale C (eps(u)-eps_free) n - gamma_s (u-g)]`.
    /// Constant modes select force components. `e_i cross (x-origin)` selects
    /// moment component i about that origin. Only the original selected
    /// embedded support is included, not additional strong box clamps.
    /// B-bar stresses retain the original physical-cell dilatation projection.
    ///
    /// Free-strain stresses are consumed from preparation, with no callback,
    /// material query, quadrature reconstruction or change to retained sources.
    /// The supplied `mode` and optional `prescribed` are pure reference laws,
    /// independent of material scales and compatible with any strong clamps.
    /// This evaluates a numerical functional and its PARTIAL derivatives;
    /// physical strain admission remains with [`Self::eigenstrain_response`].
    ///
    /// For an equilibrium sensitivity, add the adjoint's thermal load
    /// derivative from [`Self::eigenstrain_load_scale_pullback`] and any
    /// prescribed-motion lifting derivative before subtracting `z^T dK/ds u`.
    pub fn eigenstrain_embedded_reaction(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mode: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EmbeddedReaction3, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        let result = self.embedded_reaction(u, prescribed, mode, &mut checkpoint)?;
        self.subtract_eigenstress_reaction(prepared, mode, result, checkpoint)
    }

    // Accept the existing Cartesian OR already restricted adaptive derivative.
    // The correction has no displacement dependence, so no second T^T or
    // replacement support-gradient assembly belongs here.
    fn subtract_eigenstress_reaction(
        &self,
        prepared: &PreparedEigenstrain3,
        mode: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut result: EmbeddedReaction3,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EmbeddedReaction3, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        let support = self.embedded.as_ref().ok_or(ElasticityError3::Invalid(
            "thermal reaction requires original embedded support",
        ))?;
        if prepared.cells.len() != self.cells.len()
            || result.scale_gradient.len() != self.cells.len()
        {
            return Err(ElasticityError3::Invalid(
                "thermal reaction cell alignment mismatch",
            ));
        }
        for (id, (cell, free)) in self.cells.iter().zip(&prepared.cells).enumerate() {
            poll(&mut checkpoint)?;
            let selected = support.selected_points(id);
            if selected.len() != free.surface_stresses.len() {
                return Err(ElasticityError3::Invalid(
                    "thermal reaction surface alignment mismatch",
                ));
            }
            if selected.is_empty() {
                continue;
            }
            let surface = cell.rules.surface().ok_or(ElasticityError3::Invalid(
                "thermal reaction requires retained support quadrature",
            ))?;
            let mut correction = 0.0;
            for (&index, stress) in selected.iter().zip(&free.surface_stresses) {
                poll(&mut checkpoint)?;
                let point = surface
                    .points()
                    .get(index)
                    .ok_or(ElasticityError3::Invalid(
                        "thermal reaction support point is absent",
                    ))?;
                let h = mode(point.position, point.normal);
                poll(&mut checkpoint)?;
                if h.iter().chain(stress).any(|v| !v.is_finite()) {
                    return Err(ElasticityError3::Invalid(
                        "nonfinite thermal reaction mode or retained stress",
                    ));
                }
                let n = point.normal;
                let traction = [
                    stress[0] * n[0] + stress[3] * n[1] + stress[5] * n[2],
                    stress[3] * n[0] + stress[1] * n[1] + stress[4] * n[2],
                    stress[5] * n[0] + stress[4] * n[1] + stress[2] * n[2],
                ];
                correction += point.weight
                    * h.iter()
                        .zip(traction)
                        .map(|(mode, traction)| mode * traction)
                        .sum::<f64>();
                if !correction.is_finite() {
                    return Err(ElasticityError3::Invalid(
                        "thermal reaction eigenstress integral overflow",
                    ));
                }
            }
            result.scale_gradient[id] -= correction;
            result.value -= self.scales[id] * correction;
            if !result.scale_gradient[id].is_finite() || !result.value.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "thermal reaction correction overflow",
                ));
            }
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}

impl AdaptiveElasticity3 {
    /// The same thermal support reaction in independent master coordinates.
    /// The existing embedded-reaction owner reconstructs the physical field
    /// and restricts its displacement derivative once. Retained eigenstress
    /// changes only the value and direct cell-scale derivative afterwards.
    pub fn eigenstrain_embedded_reaction(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mode: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EmbeddedReaction3, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        let result = self.embedded_reaction(u, prescribed, mode, &mut checkpoint)?;
        self.raw
            .subtract_eigenstress_reaction(prepared, mode, result, checkpoint)
    }
}
