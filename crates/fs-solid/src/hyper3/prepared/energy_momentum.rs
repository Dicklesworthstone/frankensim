//! Energy/momentum-compatible forces on the original prepared P1 mesh.

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
            .energy_momentum_evaluate(before_m, after_m, false, cx)?
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
            .energy_momentum_evaluate(before_m, after_m, true, cx)?
            .1
            .expect("requested algorithmic tangent"))
    }

    fn energy_momentum_evaluate(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Option<Csr>), HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(before_m)?;
        self.problem.check_len(after_m)?;
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
            let response =
                stress::evaluate(self.problem.material.at(index), &before, &after, tangent)
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
