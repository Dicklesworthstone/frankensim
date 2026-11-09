//! Material-law controls of an already solved complete thermal adjoint.
//! This is the conduction-block contraction, not another state/adjoint solve.

use super::super::{RobinResponse, vector};
use super::{ConductionError, ConductionProblem, Cx, checked, invalid, material_is_smooth, poll};
use crate::assemble::{checked_element_conductivity, element_conduction_action};

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
    /// integrated material tensor. It does NOT insert K'(T) a second time:
    /// the full state Jacobian already accounts for it. All prescribed
    /// temperatures remain in the primal vector and their adjoint entries must
    /// be zero, retaining the material-dependent Dirichlet lift exactly once.
    /// Effective element assignments override the fallback model.
    ///
    /// Work scales with elements, vertices and crossed material knots; one
    /// bounded element-result vector, no global matrix,
    /// degree-of-freedom allocation or additional solve.
    /// max_elements bounds the output/traversal, not total allocator/RSS.
    /// Results are Estimated local model derivatives, not error certificates.
    ///
    /// # Errors
    /// Wrong field/assignment sizes, changed prescribed values, nonzero fixed
    /// adjoint entries, material kinks/domains, nonfinite arithmetic, exhausted
    /// element allowance, allocation failure, or cancellation.
    pub fn conductivity_scale_pullback(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        temperature: &[f64],
        nodal_load_adjoint: &[f64],
        max_elements: usize,
    ) -> Result<Vec<f64>, ConductionError> {
        let ne =
            material_control_elements(cx, problem, temperature, nodal_load_adjoint, max_elements)?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(ne)
            .map_err(|_| invalid("conductivity-scale result allocation refused"))?;
        for (e, tet) in problem.mesh.complex().tets.iter().enumerate() {
            poll(cx, e)?;
            let model = match problem.element_materials {
                Some(materials) => materials.model_for(e)?,
                None => problem.material,
            };
            let tensor =
                checked_element_conductivity(Some(cx), problem.mesh, e, model, temperature, false)?
                    .tensor;
            let load = element_conduction_action(problem.mesh, e, &tensor, temperature);
            let mut value = 0.0;
            for (a, &v) in tet.iter().enumerate() {
                let contribution = checked(nodal_load_adjoint[v as usize] * checked(load[a])?)?;
                value = checked(value - contribution)?;
            }
            result.push(value);
        }
        poll(cx, ne)?;
        Ok(result)
    }

    /// Pull a complete thermal adjoint back to each element's material orientation.
    ///
    /// The three entries are derivatives for right-handed infinitesimal rotations
    /// about the WORLD x, y and z axes, in objective units per radian. For
    /// R(delta) = exp([delta]_cross), the controlled law is
    /// K_e(T,delta) = R(delta) K_e(T) R(delta)^T at delta=0. The complete
    /// temperature-dependent tensor rotates; its principal conductivities,
    /// their temperature laws and positive definiteness are preserved.
    /// Sum element vectors for a region undergoing the same world-frame rotation.
    ///
    /// The integrated original tensor gives dK/d(delta_i) = Omega_i*K-K*Omega_i.
    /// This commutes with the P1 constitutive integral, including crossed table
    /// knots. The result is -lambda_e^T*(dA_e/d(delta_i))*T_e. No K'(T) is
    /// inserted again: the supplied complete nodal-load adjoint must already
    /// include material, contact, radiation and convection state feedback.
    /// Its ownership requirements are those of [Self::conductivity_scale_pullback].
    /// Isotropic laws have exactly zero rotational stiffness derivative.
    ///
    /// This rotates BULK material axes at fixed mesh geometry, reference mass,
    /// boundary/sheet/contact laws and heat sources. It does not rotate the
    /// body, modify a source claim's support, or differentiate an interface
    /// orientation law. Add explicit objective dependence separately.
    ///
    /// One 24-byte-per-element result vector is allocated. Three local
    /// contractions reuse the native matrix-free conductive action; no global
    /// matrix, extra primal solve or extra transpose solve is formed.
    ///
    /// # Errors
    /// The same field, support, smoothness, budget and cancellation refusals as
    /// [Self::conductivity_scale_pullback], plus unrepresentable derivative arithmetic.
    pub fn conductivity_rotation_pullback(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        temperature: &[f64],
        nodal_load_adjoint: &[f64],
        max_elements: usize,
    ) -> Result<Vec<[f64; 3]>, ConductionError> {
        let ne =
            material_control_elements(cx, problem, temperature, nodal_load_adjoint, max_elements)?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(ne)
            .map_err(|_| invalid("conductivity-rotation result allocation refused"))?;
        for (e, tet) in problem.mesh.complex().tets.iter().enumerate() {
            poll(cx, e)?;
            let model = match problem.element_materials {
                Some(materials) => materials.model_for(e)?,
                None => problem.material,
            };
            let tensor =
                checked_element_conductivity(Some(cx), problem.mesh, e, model, temperature, false)?
                    .tensor;
            let mut gradient = [0.0; 3];
            for (axis, derivative) in rotation_tangents(tensor).iter().enumerate() {
                let load = element_conduction_action(problem.mesh, e, derivative, temperature);
                for (a, &v) in tet.iter().enumerate() {
                    let contribution = checked(nodal_load_adjoint[v as usize] * checked(load[a])?)?;
                    gradient[axis] = checked(gradient[axis] - contribution)?;
                }
            }
            result.push(gradient);
        }
        poll(cx, ne)?;
        Ok(result)
    }
}

fn material_control_elements(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    temperature: &[f64],
    nodal_load_adjoint: &[f64],
    max_elements: usize,
) -> Result<usize, ConductionError> {
    poll(cx, 0)?;
    let n = problem.mesh.vertex_count();
    let ne = problem.mesh.element_count();
    if ne > max_elements {
        return Err(invalid(
            "conductivity pullback exceeds the element allowance",
        ));
    }
    vector(cx, temperature, n)?;
    vector(cx, nodal_load_adjoint, n)?;
    if let Some(materials) = problem.element_materials {
        materials.validate_for(problem.mesh)?;
    }
    for &(v, prescribed) in problem.boundary.dirichlet() {
        poll(cx, v)?;
        if v >= n || temperature[v] != prescribed || nodal_load_adjoint[v] != 0.0 {
            return Err(invalid(
                "conductivity pullback requires unchanged prescribed values and zero fixed-node dual",
            ));
        }
    }
    if !material_is_smooth(cx, problem, temperature)? {
        return Err(invalid(
            "conductivity pullback requires a smooth material neighborhood at the supplied field",
        ));
    }
    Ok(ne)
}

/// Exact differential of R*K*R^T at the identity, for world-axis generators.
/// No division by an eigenvalue or separation of repeated principal values:
/// isotropic tensors and repeated eigenvalues have their correct null directions.
fn rotation_tangents(k: [[f64; 3]; 3]) -> [[[f64; 3]; 3]; 3] {
    std::array::from_fn(|axis| {
        let p = (axis + 1) % 3;
        let q = (axis + 2) % 3;
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let left = if i == p {
                    -k[q][j]
                } else if i == q {
                    k[p][j]
                } else {
                    0.0
                };
                let right = if j == p {
                    -k[i][q]
                } else if j == q {
                    k[i][p]
                } else {
                    0.0
                };
                left + right
            })
        })
    })
}

#[cfg(test)]
mod tests;
