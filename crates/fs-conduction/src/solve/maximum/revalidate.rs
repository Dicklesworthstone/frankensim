//! Physical publication of a corrected field against newly marched Robin
//! references. This reuses ordinary conduction assembly and reporting; the
//! affine correction's residual is deliberately not its acceptance criterion.

use super::{
    ConductionError, ConductionProblem, ConductionSolution, Cx, DofMap, StopReason,
    ThermalInterfaces, invalid, physical, poll, refresh,
};
use crate::{ScalarField, ThermalBc, ThermalBoundary};

impl ConductionSolution {
    /// Reassemble and accept a corrected linear field against explicit uniform
    /// Robin references, returning BOTH the field/report and its boundary.
    ///
    /// Only named reference temperatures change. Face ownership, coefficients,
    /// prescribed values, sources, material assignment and matching contact
    /// stay in the supplied problem. The residual, energy balance and every
    /// Robin/contact flux are recomputed by the ordinary physical solver's
    /// report routine. Both caller-declared physical gates must pass; a small
    /// affine-model or maximum-error bound alone is insufficient.
    ///
    /// `self` supplies mesh/material metadata and historical iteration evidence,
    /// not trusted residuals. Its histories remain historical; this method does
    /// not perform a Krylov solve or invent correction iterations. The caller
    /// must retain correction work separately. The absolute watt threshold is
    /// explicit so a coupled caller can preserve its original acceptance policy.
    /// No continuum accuracy, radiation or air-law claim follows from this
    /// solid-only reassembly; the coupled caller must also check its air march.
    ///
    /// # Errors
    /// Invalid gates, mismatched fields/report, changed prescribed values,
    /// nonuniform/unknown/repeated Robin targets, nonlinear or out-of-domain
    /// material, assembly/flux refusals, failed physical gates and cancellation.
    /// Neither the original solution nor the original boundary is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn revalidate_linear_robin_temperature(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        temperature: &[f64],
        references: &[(&str, f64)],
        residual_threshold_w: f64,
        energy_relative_tolerance: f64,
    ) -> Result<(Self, ThermalBoundary), ConductionError> {
        poll(cx)?;
        if !(residual_threshold_w.is_finite() && residual_threshold_w >= 0.0)
            || !(energy_relative_tolerance.is_finite()
                && (0.0..=1.0).contains(&energy_relative_tolerance))
        {
            return Err(invalid("finite nonnegative watt threshold and energy tolerance in [0,1] required"));
        }
        let n = problem.mesh.vertex_count();
        if temperature.len() != n || self.temperature.len() != n
            || self.report.elements != problem.mesh.element_count()
        {
            return Err(invalid("corrected field and original report must match the supplied mesh"));
        }
        for (i, &value) in temperature.iter().enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            if !value.is_finite() {
                return Err(ConductionError::NonFinite { field: "corrected temperature", bits: value.to_bits() });
            }
        }
        problem.source.validate("volumetric source", n)?;
        let materials = problem.element_materials;
        if let Some(assigned) = materials { assigned.validate_for(problem.mesh)?; }
        if self.report.material_provenance != materials.map_or_else(
            || problem.material.provenance(), crate::material::ElementMaterials::provenance,
        ) || self.report.material_receipts != materials.map_or_else(
            || problem.material.receipts().len(), |assigned| assigned.receipts().len(),
        ) || self.report.element_material_identity
            != materials.map(crate::material::ElementMaterials::identity)
        {
            return Err(invalid("original material metadata does not match the supplied assignment"));
        }
        for element in 0..problem.mesh.element_count() {
            if element % 512 == 0 { poll(cx)?; }
            let model = match materials {
                Some(assigned) => assigned.model_for(element)?,
                None => problem.material,
            };
            if model.is_temperature_dependent() {
                return Err(invalid("linear field publication cannot freeze temperature-dependent conductivity"));
            }
            model.temperature_span().check(crate::assemble::element_temperature(
                problem.mesh, element, temperature,
            ))?;
        }
        if references.len() > problem.boundary.region_names().len() {
            return Err(invalid("too many Robin reference replacements"));
        }
        let mut replacements = Vec::new();
        replacements.try_reserve_exact(references.len())
            .map_err(|_| invalid("Robin reference allocation refused"))?;
        let mut seen = std::collections::BTreeSet::new();
        for &(name, reference) in references {
            poll(cx)?;
            let region = problem.boundary.region_names().iter().position(|entry| entry == name)
                .ok_or_else(|| invalid("unknown Robin reference target"))?;
            if !seen.insert(region) || !reference.is_finite() {
                return Err(invalid("Robin references must be finite and targets distinct"));
            }
            let Some(ThermalBc::Robin {
                htc: ScalarField::Uniform(htc), t_ref: ScalarField::Uniform(_),
            }) = problem.boundary.conditions().get(region) else {
                return Err(invalid("reference replacement requires an existing uniform Robin row"));
            };
            replacements.push((region, *htc, reference));
        }
        let boundary = problem.boundary.with_uniform_robin_replacements(&replacements)?;
        let dofs = DofMap::new(&boundary, n)?;
        if self.report.free_dofs != dofs.n() {
            return Err(invalid("original free-dof count differs from the supplied boundary"));
        }
        for (i, &vertex) in dofs.fixed().iter().enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            if temperature[vertex] != dofs.prescribed()[vertex] {
                return Err(invalid("corrected field changes a prescribed temperature"));
            }
        }
        let checked = physical(cx, ConductionProblem { boundary: &boundary, ..problem },
            interfaces, &dofs, temperature)?;
        if checked.final_residual > residual_threshold_w {
            return Err(ConductionError::NotConverged {
                iterations: self.report.iterations,
                residual: checked.final_residual, threshold: residual_threshold_w,
            });
        }
        if checked.energy.relative_closure() > energy_relative_tolerance {
            return Err(ConductionError::Config {
                parameter: "corrected field energy",
                what: format!("relative energy closure {} exceeds {}",
                    checked.energy.relative_closure(), energy_relative_tolerance),
            });
        }
        let mut solution = self.clone();
        solution.temperature = temperature.to_vec();
        refresh(&mut solution, checked);
        solution.report.residual_threshold = residual_threshold_w;
        solution.report.stop_reason = StopReason::ResidualTolerance;
        poll(cx)?;
        Ok((solution, boundary))
    }
}
