//! Energy/momentum-compatible forces on the original prepared P1 mesh.

use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure, algorithmic};
use fs_exec::Cx;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::{Coo, Csr};

use super::{HyperTetError, HyperTetProblem, PreparedHyperTetProblem, checkpoint};

mod stress;

struct SpatialPressure {
    before_pa: Vec<f64>,
    after_pa: Vec<f64>,
    gradient_pa_m: [f64; 3],
    clipped: bool,
}

impl<L: HyperelasticMaterial> PreparedHyperTetProblem<'_, L> {
    /// Objective discrete-gradient force between two finite configurations.
    ///
    /// Each element uses symmetric second Piola stress at the mean right
    /// Cauchy--Green metric and a strain-space work correction. Multiplication
    /// by midpoint F preserves zero internal torque at the midpoint geometry;
    /// partition of unity preserves zero internal force. The potential-work
    /// identity is subject to arithmetic and the continuous-limit threshold
    /// `||Delta C||^2 <= 64 epsilon ||Cbar||^2`. A time consumer must still gate
    /// its measured energy and nonlinear residuals. No physical damping is added.
    /// The ORIGINAL frozen material seam, including thermal references, owns
    /// every energy, stress, tangent and source-domain admission.
    pub fn energy_momentum_internal_force(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(before_m, after_m, None, false, cx)?
            .0)
    }

    /// Analytic derivative of the discrete-gradient force with respect to the
    /// AFTER displacement, with the initial configuration held fixed.
    /// A solver whose unknown is midpoint displacement must multiply this
    /// action by two. The tangent is generally nonsymmetric.
    pub fn energy_momentum_tangent(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(before_m, after_m, None, true, cx)?
            .1
            .expect("requested algorithmic tangent"))
    }

    /// Material internal force MINUS uniform closed-boundary pressure force,
    /// using one objective discrete gradient of W(F)+p det(F). Positive signed
    /// pressure acts inward. Surface identity/closure remain owned by the native
    /// pressure surface; its endpoint volume work must still be audited by the
    /// time consumer. This finite-step algorithmic volume potential need not
    /// equal pointwise midpoint face traction or cancel at interior nodes;
    /// it approaches the native boundary load as the endpoints coalesce.
    /// Per-face nonconservative pressure is outside this entry.
    pub fn energy_momentum_force_with_uniform_pressure(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        surface: &ClosedTetPressureSurface,
        pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(
                before_m,
                after_m,
                Some((surface, TetPressure::Uniform(pressure_pa))),
                false,
                cx,
            )?
            .0)
    }

    /// Endpoint derivative of the combined material/pressure internal force.
    /// As with the unpressurized tangent, midpoint unknowns require a factor two.
    pub fn energy_momentum_tangent_with_uniform_pressure(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        surface: &ClosedTetPressureSurface,
        pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(
                before_m,
                after_m,
                Some((surface, TetPressure::Uniform(pressure_pa))),
                true,
                cx,
            )?
            .1
            .expect("requested algorithmic tangent"))
    }

    /// Material force minus conservative closed-surface pressure, including
    /// affine Eulerian fields and prescribed planar free surfaces.
    ///
    /// Each tetrahedron uses its exact native potential V*m(p_nodes), with
    /// discrete gradient m_bar*DG(V)+V_bar*DG(m). The volume term retains the
    /// original objective strain-space method; field forces are parallel to
    /// the spatial pressure gradient. Thus unforced transverse linear momentum
    /// and angular momentum about that gradient are preserved at midpoint,
    /// while pressure supplies buoyancy along the gradient axis.
    /// The time consumer must still gate native endpoint work and residuals.
    /// Uniform fields retain the existing uniform algorithm exactly. A held
    /// PerFace field has no assumed potential and is refused by this entry.
    pub fn energy_momentum_force_with_pressure<'p>(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(
                before_m,
                after_m,
                Some((surface, pressure.into())),
                false,
                cx,
            )?
            .0)
    }

    /// Derivative of the conservative algorithmic force with respect to AFTER
    /// displacement. Differentiates wet-volume moments and coordinate secants
    /// analytically using the native scalar AD; midpoint unknowns need factor 2.
    pub fn energy_momentum_tangent_with_pressure<'p>(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        Ok(self
            .energy_momentum_evaluate(
                before_m,
                after_m,
                Some((surface, pressure.into())),
                true,
                cx,
            )?
            .1
            .expect("requested algorithmic tangent"))
    }

    fn energy_momentum_evaluate(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Option<Csr>), HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(before_m)?;
        self.problem.check_len(after_m)?;
        let mut spatial = None;
        let pressure_pa = if let Some((surface, value)) = pressure {
            if !surface.matches_reference(self.problem.nodes_m, self.problem.tetrahedra) {
                return Err(HyperTetError::InvalidInput {
                    what: "algorithmic pressure surface does not match the material reference mesh"
                        .into(),
                });
            }
            surface.validate_pressure(value, cx)?;
            match value {
                TetPressure::Uniform(p) => p,
                TetPressure::PerFace(_) => return Err(HyperTetError::InvalidInput {
                    what: "energy-momentum pressure requires a conservative uniform, affine or free-surface field".into(),
                }),
                TetPressure::Affine { origin_m, pressure_pa, gradient_pa_m }
                | TetPressure::FreeSurface { origin_m, pressure_pa, gradient_pa_m } => {
                    let clipped = matches!(value, TetPressure::FreeSurface { .. });
                    if gradient_pa_m == [0.0; 3] {
                        if clipped { pressure_pa.max(0.0) } else { pressure_pa }
                    } else {
                        let mut before_pa = Vec::with_capacity(self.problem.nodes_m.len());
                        let mut after_pa = Vec::with_capacity(self.problem.nodes_m.len());
                        for node in 0..self.problem.nodes_m.len() {
                            checkpoint(cx)?;
                            before_pa.push(surface.spatial_pressure_at(before_m, node, origin_m, pressure_pa, gradient_pa_m)?);
                            after_pa.push(surface.spatial_pressure_at(after_m, node, origin_m, pressure_pa, gradient_pa_m)?);
                        }
                        spatial = Some(SpatialPressure { before_pa, after_pa, gradient_pa_m, clipped });
                        0.0
                    }
                }
            }
        } else {
            0.0
        };
        let mut force = vec![0.0; after_m.len()];
        let mut matrix = tangent.then(|| Coo::new(after_m.len(), after_m.len()));
        for (index, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let before = HyperTetProblem::<L>::deformation_gradient(element, before_m);
            let after = HyperTetProblem::<L>::deformation_gradient(element, after_m);
            let mut jacobians = [0.0; 2];
            for (endpoint, f) in [&before, &after].into_iter().enumerate() {
                jacobians[endpoint] =
                    HyperTetProblem::<L>::deformation_jacobian(f).map_err(|what| {
                        HyperTetError::MaterialRefused {
                            element: index,
                            what,
                        }
                    })?;
            }
            let mean = spatial
                .as_ref()
                .map(|field| {
                    algorithmic::evaluate(
                        element.nodes.map(|node| field.before_pa[node]),
                        element.nodes.map(|node| field.after_pa[node]),
                        field.clipped,
                        tangent,
                    )
                })
                .transpose()?;
            let held_mean = mean.as_ref().map_or(pressure_pa, |mean| {
                0.5 * mean.before_pa + 0.5 * mean.after_pa
            });
            let volume_mean = element.volume * (0.5 * jacobians[0] + 0.5 * jacobians[1]);
            let response = stress::evaluate_with_uniform_pressure(
                self.problem.material.at(index),
                &before,
                &after,
                held_mean,
                tangent,
            )
            .map_err(|error| HyperTetError::MaterialRefused {
                element: index,
                what: error.to_string(),
            })?;
            for (a, &node) in element.nodes.iter().enumerate() {
                for i in 0..3 {
                    let value = (0..3).fold(0.0, |s, j| {
                        response.piola[3 * i + j].mul_add(element.gradients[a][j], s)
                    });
                    force[3 * node + i] = element.volume.mul_add(value, force[3 * node + i]);
                    if let (Some(field), Some(mean)) = (&spatial, &mean) {
                        force[3 * node + i] +=
                            volume_mean * mean.weights[a] * field.gradient_pa_m[i];
                    }
                }
            }
            if let (Some(matrix), Some(tangent)) = (&mut matrix, response.tangent_after) {
                let volume_piola = if mean
                    .as_ref()
                    .is_some_and(|m| m.after_derivative.iter().any(|d| *d != 0.0))
                {
                    stress::volume_stress(self.problem.material.at(index), &before, &after)
                        .map_err(|error| HyperTetError::MaterialRefused {
                            element: index,
                            what: error.to_string(),
                        })?
                } else {
                    [0.0; 9]
                };
                let after_cofactor = if mean
                    .as_ref()
                    .is_some_and(|m| m.weights.iter().any(|w| *w != 0.0))
                {
                    stress::cofactor(&after)
                } else {
                    [0.0; 9]
                };
                let volume_forces: [[f64; 3]; 4] = core::array::from_fn(|a| {
                    core::array::from_fn(|i| {
                        (0..3).fold(0.0, |s, j| {
                            volume_piola[3 * i + j].mul_add(element.gradients[a][j], s)
                        })
                    })
                });
                let after_volume_derivatives: [[f64; 3]; 4] = core::array::from_fn(|b| {
                    core::array::from_fn(|k| {
                        (0..3).fold(0.0, |s, l| {
                            after_cofactor[3 * k + l].mul_add(element.gradients[b][l], s)
                        })
                    })
                });
                for a in 0..4 {
                    for i in 0..3 {
                        for b in 0..4 {
                            for k in 0..3 {
                                let mut value = 0.0;
                                for j in 0..3 {
                                    for l in 0..3 {
                                        value = (tangent[3 * i + j][3 * k + l]
                                            * element.gradients[a][j])
                                            .mul_add(element.gradients[b][l], value);
                                    }
                                }
                                let mut entry = element.volume * value;
                                if let (Some(field), Some(mean)) = (&spatial, &mean) {
                                    entry += element.volume
                                        * volume_forces[a][i]
                                        * (0.5 * mean.after_derivative[b] * field.gradient_pa_m[k])
                                        + 0.5
                                            * element.volume
                                            * after_volume_derivatives[b][k]
                                            * mean.weights[a]
                                            * field.gradient_pa_m[i]
                                        + volume_mean
                                            * (mean.tangent_after[a][b] * field.gradient_pa_m[k])
                                            * field.gradient_pa_m[i];
                                }
                                matrix.push(
                                    3 * element.nodes[a] + i,
                                    3 * element.nodes[b] + k,
                                    entry,
                                );
                            }
                        }
                    }
                }
            }
        }
        let matrix = matrix.map(|coo| coo.assemble());
        if force.iter().any(|v| !v.is_finite())
            || matrix.as_ref().is_some_and(|m| {
                (0..m.nrows()).any(|row| m.row(row).1.iter().any(|v| !v.is_finite()))
            })
        {
            return Err(HyperTetError::InvalidInput {
                what: "unrepresentable discrete-gradient force or tangent".into(),
            });
        }
        checkpoint(cx)?;
        Ok((force, matrix))
    }
}
