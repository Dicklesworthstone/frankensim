//! Material-law controls of an already solved complete thermal adjoint.
//! This is the conduction-block contraction, not another state/adjoint solve.

use super::{ConductionError, ConductionProblem, Cx, checked, invalid,
    material_is_smooth, poll};
use super::super::{RobinResponse, vector};
use crate::assemble::{element_stiffness, element_temperature};

impl RobinResponse {
    /// Pull a complete nodal-load adjoint back to independent conductivity scales.
    ///
    /// In element e the control is `K_e(T) -> exp(p_e) K_e(T)`, evaluated at
    /// p_e=0. Return dJ/dp_e in element order. This also equals dJ/ds_e at s_e=1
    /// for a multiplicative scale s_e of the retained material law. Summing
    /// entries over a region gives its shared material-law scale derivative.
    /// The WHOLE tensor/curve scales; this is not a derivative with respect to
    /// an absolute scalar conductivity, an individual tensor entry, a table
    /// knot, material selection, heat capacity, geometry or contact resistance.
    ///
    /// The caller supplies the unchanged physical field and the SAME total
    /// load adjoint produced by its complete state solve. That dual must
    /// already contain smooth k(T), contact, radiation, natural convection and
    /// air feedback as applicable. This function cannot attest that external
    /// binding or its full residual; a frozen-boundary dual is not made total
    /// by this contraction. Objectives with an explicit material dependence
    /// must add that direct term separately.
    ///
    /// The fixed-state partial uses the production element stiffness at the
    /// production element temperature. It does NOT insert K'(T) a second time:
    /// the full state Jacobian already accounts for it. All prescribed
    /// temperatures remain in the primal vector and their adjoint entries must
    /// be zero, retaining the material-dependent Dirichlet lift exactly once.
    /// Effective element assignments override the fallback model.
    ///
    /// O(elements + vertices) work, one bounded element-result vector; no global matrix,
    /// degree-of-freedom allocation or additional solve.
    /// max_elements bounds the output/traversal, not total allocator/RSS.
    /// Results are Estimated local model derivatives, not error certificates.
    ///
    /// # Errors
    /// Wrong field/assignment sizes, changed prescribed values, nonzero fixed
    /// adjoint entries, material kinks/domains, nonfinite arithmetic, exhausted
    /// element allowance, allocation failure, or cancellation.
    pub fn conductivity_scale_pullback(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
        nodal_load_adjoint: &[f64], max_elements: usize,
    ) -> Result<Vec<f64>, ConductionError> {
        poll(cx, 0)?;
        let n = problem.mesh.vertex_count();
        let ne = problem.mesh.element_count();
        if ne > max_elements {
            return Err(invalid("conductivity-scale pullback exceeds the element allowance"));
        }
        vector(cx, temperature, n)?;
        vector(cx, nodal_load_adjoint, n)?;
        if let Some(materials) = problem.element_materials {
            materials.validate_for(problem.mesh)?;
        }
        for &(v, prescribed) in problem.boundary.dirichlet() {
            poll(cx, v)?;
            if v >= n || temperature[v] != prescribed || nodal_load_adjoint[v] != 0.0 {
                return Err(invalid("conductivity-scale pullback requires the unchanged prescribed values and zero fixed-node dual"));
            }
        }
        if !material_is_smooth(cx, problem, temperature)? {
            return Err(invalid("conductivity-scale pullback cannot choose a total derivative at a material kink or validity endpoint"));
        }
        let mut result = Vec::new();
        result.try_reserve_exact(ne)
            .map_err(|_| invalid("conductivity-scale result allocation refused"))?;
        for (e, tet) in problem.mesh.complex().tets.iter().enumerate() {
            poll(cx, e)?;
            let model = match problem.element_materials {
                Some(materials) => materials.model_for(e)?,
                None => problem.material,
            };
            let sample = checked(element_temperature(problem.mesh, e, temperature))?;
            let tensor = model.tensor_at(sample)?;
            let stiffness = element_stiffness(problem.mesh, e, &tensor);
            let mut value = 0.0;
            for (a, &v) in tet.iter().enumerate() {
                for (b, &w) in tet.iter().enumerate() {
                    let contribution = checked(checked(nodal_load_adjoint[v as usize]
                        * checked(stiffness[a][b])?)? * temperature[w as usize])?;
                    value = checked(value - contribution)?;
                }
            }
            result.push(value);
        }
        poll(cx, ne)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
