//! Energy/momentum-compatible forces on the original prepared P1 mesh.

use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_exec::Cx;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::{Coo, Csr};

use super::{HyperTetError, HyperTetProblem, PreparedHyperTetProblem, checkpoint};

mod stress;

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
            .energy_momentum_evaluate(before_m, after_m, Some((surface, pressure_pa)), false, cx)?
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
            .energy_momentum_evaluate(before_m, after_m, Some((surface, pressure_pa)), true, cx)?
            .1
            .expect("requested algorithmic tangent"))
    }

    fn energy_momentum_evaluate(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        pressure: Option<(&ClosedTetPressureSurface, f64)>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Option<Csr>), HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(before_m)?;
        self.problem.check_len(after_m)?;
        let pressure_pa = if let Some((surface, value)) = pressure {
            if !surface.matches_reference(self.problem.nodes_m, self.problem.tetrahedra) {
                return Err(HyperTetError::InvalidInput {
                    what: "algorithmic pressure surface does not match the material reference mesh"
                        .into(),
                });
            }
            surface.validate_pressure(TetPressure::Uniform(value), cx)?;
            value
        } else {
            0.0
        };
        let mut force = vec![0.0; after_m.len()];
        let mut matrix = tangent.then(|| Coo::new(after_m.len(), after_m.len()));
        for (index, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let before = HyperTetProblem::<L>::deformation_gradient(element, before_m);
            let after = HyperTetProblem::<L>::deformation_gradient(element, after_m);
            for f in [&before, &after] {
                HyperTetProblem::<L>::deformation_jacobian(f).map_err(|what| {
                    HyperTetError::MaterialRefused {
                        element: index,
                        what,
                    }
                })?;
            }
            let response = stress::evaluate_with_uniform_pressure(
                self.problem.material.at(index),
                &before,
                &after,
                pressure_pa,
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
                }
            }
            if let (Some(matrix), Some(tangent)) = (&mut matrix, response.tangent_after) {
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
                                matrix.push(
                                    3 * element.nodes[a] + i,
                                    3 * element.nodes[b] + k,
                                    element.volume * value,
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
