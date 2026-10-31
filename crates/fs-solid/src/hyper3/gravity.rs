//! Conservative gravity sampled on the CURRENT P1 material configuration.
//!
//! The four equal-weight tetrahedral points integrate quadratics exactly.
//! With fixed reference mass weights `w`, the SAME discrete potential gives
//! `U = sum w phi(x_q)`, `f_ai = sum w N_a g_i(x_q)` and
//! `df_ai/du_bj = sum w N_a N_b dg_i/dx_j`. Thus uniform and affine gravity
//! reproduce the native consistent mass projection, including supported DOFs.
//! General smooth fields require spatial refinement; a four-point rule does
//! not resolve a nearby singular source or certify a quadrature error bound.
//!
//! Gravity is an external conservative body field, not a current-area follower
//! pressure and not a force density multiplied by the changing material volume.
//! The callback must be deterministic and stationary, and must supply mutually
//! consistent potential, acceleration and gradient in the mesh's world axes.

use fs_exec::Cx;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::{Coo, Csr};

use super::{Element, HyperTetError, HyperTetProblem, checkpoint};
use crate::linear3::p1_tet_mass_scale;

/// One conservative field evaluation in fixed world coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GravitySample3 {
    /// Specific potential [J/kg], with one fixed datum for the entire field.
    pub potential_j_kg: f64,
    /// Acceleration [m/s²], equal to minus the gradient of `potential_j_kg`.
    pub acceleration_m_s2: [f64; 3],
    /// Spatial acceleration derivative [1/s²]: entry `(i,j)` is `dg_i/dx_j`.
    /// It must be symmetric to roundoff for this conservative field interface.
    pub acceleration_gradient_s2: [[f64; 3]; 3],
}

/// A borrowed stationary field at a current world position [m].
///
/// Refuse outside the original field's domain; do not clamp a singularity or
/// extrapolate source data. Derivative consistency is the field owner's
/// obligation, not a property established by the Rust callback type.
pub type GravityField3<'a> = dyn Fn([f64; 3]) -> Result<GravitySample3, HyperTetError> + 'a;

/// Gravity integrated with the original tetrahedral reference mass inventory.
#[derive(Debug)]
pub struct TetGravityLoad3 {
    /// Full node-major world forces [N], including prescribed coordinates.
    pub force_n: Vec<f64>,
    /// Derivative of the EXTERNAL force [N/m]; subtract from internal stiffness.
    pub force_tangent_n_m: Csr,
    /// Gravitational potential [J], separate from elastic strain energy.
    pub potential_energy_j: f64,
}

pub(super) struct GravityEvaluation {
    pub force_n: Vec<f64>,
    pub tangent: Option<Csr>,
    pub potential_energy_j: f64,
}

fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput { what: what.into() }
}

fn finite(value: f64) -> Result<f64, HyperTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("gravity sample or integral is not representable"))
    }
}

impl GravitySample3 {
    fn validate(self) -> Result<Self, HyperTetError> {
        finite(self.potential_j_kg)?;
        for value in self.acceleration_m_s2 {
            finite(value)?;
        }
        let mut scale = f64::MIN_POSITIVE;
        for row in self.acceleration_gradient_s2 {
            for value in row {
                finite(value)?;
                scale = scale.max(value.abs());
            }
        }
        for i in 0..3 {
            for j in 0..i {
                if (self.acceleration_gradient_s2[i][j] - self.acceleration_gradient_s2[j][i]).abs()
                    > 64.0 * f64::EPSILON * scale
                {
                    return Err(invalid(
                        "conservative gravity needs a symmetric acceleration gradient",
                    ));
                }
            }
        }
        Ok(self)
    }
}

/// One quarter of each element's fixed reference mass, shared by all four
/// integration points. The scalar is the native consistent P1 mass kernel.
pub(super) fn reference_weights(
    elements: &[Element],
    reference_density_kg_m3: &[f64],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, HyperTetError> {
    checkpoint(cx)?;
    if reference_density_kg_m3.len() != elements.len() {
        return Err(invalid(
            "gravity needs one reference density per tetrahedron",
        ));
    }
    let mut weights = Vec::with_capacity(elements.len());
    for (element, &density) in elements.iter().zip(reference_density_kg_m3) {
        checkpoint(cx)?;
        if !density.is_finite() || density <= 0.0 {
            return Err(invalid("gravity needs finite positive reference densities"));
        }
        let weight = finite(5.0 * p1_tet_mass_scale(element.volume, density))?;
        if weight <= 0.0 {
            return Err(invalid("gravity reference mass is unrepresentable"));
        }
        weights.push(weight);
    }
    Ok(weights)
}

// N_i values of the symmetric degree-two tetrahedral rule. Each point's
// weight is V/4: sum w N_i = V/4 and sum w N_i N_j = V(1+delta_ij)/20.
const LARGE: f64 = 0.585_410_196_624_968_5;
const SMALL: f64 = 0.138_196_601_125_010_5;

#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub(super) fn evaluate(
    nodes_m: &[[f64; 3]],
    elements: &[Element],
    mass_weights: &[f64],
    displacement_m: &[f64],
    field: &GravityField3<'_>,
    with_tangent: bool,
    cx: &Cx<'_>,
) -> Result<GravityEvaluation, HyperTetError> {
    checkpoint(cx)?;
    let dofs = 3 * nodes_m.len();
    if displacement_m.len() != dofs || displacement_m.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "gravity needs one finite xyz displacement per reference node",
        ));
    }
    let mut force_n = vec![0.0; dofs];
    let mut potential_energy_j = 0.0;
    let mut entries = with_tangent.then(|| Coo::new(dofs, dofs));
    for (element, &weight) in elements.iter().zip(mass_weights) {
        checkpoint(cx)?;
        let mut local = [[0.0; 12]; 12];
        let origin = element.nodes[0];
        for point in 0..4 {
            let shape = std::array::from_fn::<_, 4, _>(|a| if a == point { LARGE } else { SMALL });
            // Relative interpolation preserves common translation in the
            // geometry instead of adding four large absolute coordinates.
            let mut position = [0.0; 3];
            for i in 0..3 {
                position[i] = finite(nodes_m[origin][i] + displacement_m[3 * origin + i])?;
                for a in 1..4 {
                    let node = element.nodes[a];
                    let relative = finite(
                        (nodes_m[node][i] - nodes_m[origin][i])
                            + (displacement_m[3 * node + i] - displacement_m[3 * origin + i]),
                    )?;
                    position[i] = finite(shape[a].mul_add(relative, position[i]))?;
                }
            }
            let sample = field(position)?.validate()?;
            potential_energy_j = finite(weight.mul_add(sample.potential_j_kg, potential_energy_j))?;
            for a in 0..4 {
                let share = weight * shape[a];
                for i in 0..3 {
                    let dof = 3 * element.nodes[a] + i;
                    force_n[dof] =
                        finite(share.mul_add(sample.acceleration_m_s2[i], force_n[dof]))?;
                }
                if with_tangent {
                    for b in 0..4 {
                        let product = share * shape[b];
                        for i in 0..3 {
                            for j in 0..3 {
                                local[3 * a + i][3 * b + j] = finite(product.mul_add(
                                    sample.acceleration_gradient_s2[i][j],
                                    local[3 * a + i][3 * b + j],
                                ))?;
                            }
                        }
                    }
                }
            }
        }
        if let Some(entries) = &mut entries {
            for a in 0..4 {
                for b in 0..4 {
                    for i in 0..3 {
                        for j in 0..3 {
                            let value = local[3 * a + i][3 * b + j];
                            if value != 0.0 {
                                entries.push(
                                    3 * element.nodes[a] + i,
                                    3 * element.nodes[b] + j,
                                    value,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    let tangent = entries.map(|entries| entries.assemble());
    if let Some(matrix) = &tangent {
        for row in 0..dofs {
            if row % 4096 == 0 {
                checkpoint(cx)?;
            }
            for &value in matrix.row(row).1 {
                finite(value)?;
            }
        }
    }
    checkpoint(cx)?;
    Ok(GravityEvaluation {
        force_n,
        tangent,
        potential_energy_j,
    })
}

impl<L: HyperelasticMaterial> HyperTetProblem<'_, L> {
    /// Sample current-configuration gravity with the native reference mass.
    ///
    /// Densities are exactly the element-ordered reference densities used by
    /// [`Self::reference_mass_matrix`]; they must remain fixed during motion.
    /// The field has no implicit Earth preset. This operation assembles a load;
    /// native equilibrium separately admits the current constitutive state.
    ///
    /// # Errors
    /// Invalid reference problem, density/displacement/field, overflow or cancellation.
    pub fn gravity_load(
        &self,
        reference_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetGravityLoad3, HyperTetError> {
        self.prepare(cx)?
            .gravity_load(reference_density_kg_m3, displacement_m, field, cx)
    }
}
