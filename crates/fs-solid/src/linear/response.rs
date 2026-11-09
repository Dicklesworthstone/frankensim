//! Observe the original P1/Q1 linear-elastic weak form before support elimination.
//!
//! Stress and stored energy use the same Standard or B-bar strain columns as
//! the stiffness. B-bar values are explicitly assumed strains; the compatible
//! in-plane strain remains available separately. This owner projects the TWO
//! in-plane normal strains, preserving epsilon_zz=0 in plane strain. It is not
//! the distinct full 3-D projection used by the generic plasticity frontend.
//!
//! All forces and integrated energies are per reference thickness, matching
//! `LinearProblem`: multiply by a physical thickness in metres to obtain N and
//! J. Quadrature strain and stress are dimensionless and Pa. Recovery does not
//! solve equilibrium, certify a material strain domain, or validate small
//! rotations; the caller must choose appropriate linear-elastic kinematics.

use super::{LinearProblem, PlaneKind, SolidError, strain_column};
use crate::mesh2::{quad_points, shapes_at};
use fs_exec::Cx;
use fs_material::Voigt;

/// One original integration point, without nodal extrapolation or averaging.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearElasticPoint {
    /// Index in the original mesh connectivity.
    pub element: usize,
    /// Local index in the original `mesh2::quad_points` order.
    pub local_point: usize,
    /// Original physical position of this quadrature sample.
    pub reference_position_m: [f64; 2],
    /// Positive reference-area weight, including the Jacobian.
    pub area_weight_m2: f64,
    /// Compatible `[epsilon_xx, epsilon_yy, epsilon_xy]` (TENSOR shear).
    pub compatible_in_plane_strain: [f64; 3],
    /// Constitutive `[xx,yy,zz,xy,yz,zx]` with tensor shear. Plane stress
    /// retains its relaxed epsilon_zz; B-bar retains the projected xx/yy.
    pub strain: Voigt,
    /// Work-conjugate Cauchy stress [Pa] in the same tensor ordering.
    /// Plane strain retains sigma_zz, while plane stress has sigma_zz=0.
    pub stress_pa: Voigt,
    /// Original constitutive storage, `sigma:epsilon/2`.
    pub stored_energy_density_j_m3: f64,
}

/// Full physical response of a supplied displacement; no identity support rows.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearElasticResponse {
    /// Element-major samples in the stiffness assembly's original order.
    pub quadrature: Vec<LinearElasticPoint>,
    /// Reference-area integral of original elastic energy density.
    pub stored_energy_j_m: f64,
    /// Internal minus applied force [N/m], node-major xy. Free rows retain
    /// their actual equilibrium error instead of being silently zeroed.
    pub nodal_residual_n_m: Vec<[f64; 2]>,
    /// `(2*node+component, N/m)` on exactly the original prescribed/symmetry
    /// rows. These are support reactions when the supplied field is an
    /// equilibrium satisfying the original boundary data; otherwise they are
    /// the required additional forces at those rows. Applied body and traction
    /// loads on supports are included before taking this difference.
    pub support_reactions_n_m: Vec<(usize, f64)>,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

impl LinearProblem<'_> {
    /// Recover original quadrature strain, stress, storage and support forces.
    ///
    /// Does not mutate the problem or displacement. Cancellation is polled
    /// between elements and loaded edges; a refusal publishes no partial result.
    /// Material, geometry and support admission are shared with `solve`.
    ///
    /// # Errors
    /// Invalid mesh/material/supports, incomplete or nonfinite displacement,
    /// nonfinite applied loads, unrepresentable response, or cancellation.
    pub fn response(
        &self,
        displacement_m: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<LinearElasticResponse, SolidError> {
        self.response_scaled(displacement_m, None, cx)
    }

    /// Recover the same response with the exact per-element multipliers used
    /// by `solve_with_element_scales`. Loads and supports are not scaled.
    ///
    /// # Errors
    /// The same refusals as `response`, plus invalid element multipliers.
    pub fn response_with_element_scales(
        &self,
        displacement_m: &[[f64; 2]],
        scales: &[f64],
        cx: &Cx<'_>,
    ) -> Result<LinearElasticResponse, SolidError> {
        self.response_scaled(displacement_m, Some(scales), cx)
    }

    fn response_scaled(
        &self,
        displacement: &[[f64; 2]],
        scales: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<LinearElasticResponse, SolidError> {
        let checkpoint = || cx.checkpoint().map_err(|_| SolidError::Cancelled);
        checkpoint()?;
        let (lambda, mu) = self.validate(scales)?;
        if displacement.len() != self.mesh.nodes.len()
            || displacement.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite xy displacement per original node is required",
            ));
        }
        let fixed = self.fixed_values()?;
        let mut quadrature = Vec::new();
        let mut energy = 0.0;
        let mut residual = vec![[0.0; 2]; displacement.len()];
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            checkpoint()?;
            let scale = scales.map_or(1.0, |values| values[element]);
            let (lambda, mu) = (lambda * scale, mu * scale);
            let mean = self.mean_gradients(conn);
            for (local_point, (xi, eta, w)) in quad_points(conn.len()).into_iter().enumerate() {
                let (shape, gradients, determinant) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let area = w * determinant;
                let mut compatible = [0.0; 3];
                let mut assumed = [0.0; 3];
                let mut position = [0.0; 2];
                for (a, &node) in conn.iter().enumerate() {
                    for component in 0..2 {
                        // Partition of unity makes any uniform translation
                        // strain-free, including on skew elements.
                        let du = displacement[node][component] - displacement[conn[0]][component];
                        let b = strain_column(&gradients, None, a, component);
                        let bbar = strain_column(&gradients, mean.as_deref(), a, component);
                        for i in 0..3 {
                            compatible[i] += b[i] * du;
                            assumed[i] += bbar[i] * du;
                        }
                        position[component] += shape[a] * self.mesh.nodes[node][component];
                    }
                }
                let trace = assumed[0] + assumed[1];
                let difference = assumed[0] - assumed[1];
                let mean_stress = (lambda + mu) * trace;
                let deviator = mu * difference;
                let sxx = mean_stress + deviator;
                let syy = mean_stress - deviator;
                let sxy = mu * assumed[2];
                let (ezz, szz) = match self.plane {
                    PlaneKind::Strain => (0.0, lambda * trace),
                    PlaneKind::Stress => (-self.poisson / (1.0 - self.poisson) * trace, 0.0),
                };
                // Positive modal split of 1/2 sigma:epsilon, without subtracting
                // large normal-stress products during nearly isochoric loading.
                let density =
                    0.5 * (mean_stress * trace + deviator * difference + sxy * assumed[2]);
                let strain = [assumed[0], assumed[1], ezz, 0.5 * assumed[2], 0.0, 0.0];
                let stress = [sxx, syy, szz, sxy, 0.0, 0.0];
                let body = self.point_body_force(conn, &shape).unwrap_or([0.0; 2]);
                if !(area.is_finite() && area > 0.0 && density.is_finite())
                    || strain
                        .iter()
                        .chain(&stress)
                        .chain(&body)
                        .chain(&compatible)
                        .any(|v| !v.is_finite())
                {
                    return Err(invalid(
                        "unrepresentable quadrature strain, stress, storage or body force",
                    ));
                }
                energy += area * density;
                for (a, &node) in conn.iter().enumerate() {
                    for component in 0..2 {
                        let b = strain_column(&gradients, mean.as_deref(), a, component);
                        let internal = b[0] * sxx + b[1] * syy + b[2] * sxy;
                        residual[node][component] += area * (internal - shape[a] * body[component]);
                    }
                }
                quadrature.push(LinearElasticPoint {
                    element,
                    local_point,
                    reference_position_m: position,
                    area_weight_m2: area,
                    compatible_in_plane_strain: [compatible[0], compatible[1], 0.5 * compatible[2]],
                    strain,
                    stress_pa: stress,
                    stored_energy_density_j_m3: density,
                });
            }
        }
        self.traction_forces(&checkpoint, |dof, force| {
            residual[dof / 2][dof % 2] -= force
        })?;
        if !energy.is_finite() || residual.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable integrated elastic response"));
        }
        let support_reactions_n_m = fixed
            .keys()
            .map(|&dof| (dof, residual[dof / 2][dof % 2]))
            .collect();
        checkpoint()?;
        Ok(LinearElasticResponse {
            quadrature,
            stored_energy_j_m: energy,
            nodal_residual_n_m: residual,
            support_reactions_n_m,
        })
    }
}
