//! Exact design transposes of the original prepared free-strain load and stress.
//!
//! A stiffness scale multiplies both K and the equivalent eigenstrain load.
//! For K(s) u = f + b_free(s), the equilibrium contribution to an objective
//! gradient is z^T db_free/ds - z^T dK/ds u, where K^T z = dJ/du.
//! The objective's explicit scale derivative is separate. In particular,
//! compliance with the COMPLETE scale-dependent load has TWO load-work terms;
//! the fixed-load shortcut -u^T dK/ds u is not its derivative.

use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::PreparedEigenstrain3;
use crate::elastic3::stress::BulkStressPullback3;
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3};

fn checkpoint(poll: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if poll().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

impl CutElasticity3 {
    /// Exact `dual^T db_free/dscale_c` for each original active cell.
    ///
    /// The prepared column already contains the bulk eigenstress term and the
    /// negative embedded Nitsche eigenstress traction. No stiffness scale is
    /// applied a second time. Homogeneous strong clamps contribute no work.
    /// The dual is algebraic: a physical strain-domain limit does not apply.
    /// Geometry, material tensors, free strain and boundary assignment are held
    /// fixed. Combine with [`Self::scale_bilinear_forms`] and the chosen
    /// objective's direct derivative for an equilibrium sensitivity.
    pub fn eigenstrain_load_scale_pullback(
        &self,
        prepared: &PreparedEigenstrain3,
        dual: &[f64],
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.eigenstrain_load_scale_pullback_mapped(
            prepared,
            dual,
            self.n(),
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            poll,
        )
    }

    pub(in crate::elastic3) fn eigenstrain_load_scale_pullback_mapped(
        &self,
        prepared: &PreparedEigenstrain3,
        dual: &[f64],
        n: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        checkpoint(&mut poll)?;
        self.validate_eigenstrain(prepared)?;
        if dual.len() != n {
            return Err(ElasticityError3::Invalid("eigenstrain dual shape mismatch"));
        }
        for chunk in dual.chunks(192) {
            checkpoint(&mut poll)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid("nonfinite eigenstrain dual"));
            }
        }
        let mut work = Vec::with_capacity(self.cells.len());
        for (cell, source) in self.cells.iter().zip(&prepared.cells) {
            checkpoint(&mut poll)?;
            let dual = local(cell, dual);
            if dual.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain reconstructed dual overflow",
                ));
            }
            let value: f64 = dual.iter().zip(source.force).map(|(z, b)| z * b).sum();
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain load pullback overflow",
                ));
            }
            work.push(value);
        }
        checkpoint(&mut poll)?;
        Ok(work)
    }

    /// Exact partial derivatives of a scalar objective of thermal stress.
    ///
    /// One cotangent `[xx, yy, zz, xy, yz, zx]` is required per ORIGINAL bulk
    /// point. These pair with physical stress components directly; callers
    /// include any quadrature weights in their cotangents. `displacement`
    /// contains `dJ/du`, and `scales` contains `dJ/dscale` at fixed u and fixed
    /// free strain. The latter contracts `C (eps(u) - eps_free)`, not `C eps(u)`.
    /// Both strain projection and primal-domain admission use the original
    /// eigenstrain response owner. Directions never consume a material strain
    /// allowance. This is not a geometry or source-temperature derivative.
    pub fn eigenstrain_stress_pullback(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        derivatives: &[[f64; 6]],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        checkpoint(&mut poll)?;
        for derivative in derivatives {
            checkpoint(&mut poll)?;
            if derivative.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "nonfinite eigenstrain stress cotangent",
                ));
            }
        }
        let points = self.eigenstrain_response(
            prepared,
            u,
            maximum_kinematic_strain,
            derivatives.len(),
            &mut poll,
        )?;
        if points.len() != derivatives.len() {
            return Err(ElasticityError3::Invalid(
                "one eigenstrain stress cotangent per bulk point required",
            ));
        }
        let mut scales = vec![0.0; self.cells.len()];
        let mut scaled = Vec::with_capacity(derivatives.len());
        for (point, derivative) in points.iter().zip(derivatives) {
            checkpoint(&mut poll)?;
            scales[point.cell] += derivative
                .iter()
                .zip(point.reference_stress)
                .map(|(d, stress)| d * stress)
                .sum::<f64>();
            let cotangent = derivative.map(|value| self.scales[point.cell] * value);
            if !scales[point.cell].is_finite() || cotangent.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain stress pullback overflow",
                ));
            }
            scaled.push(cotangent);
        }
        // Subtracting prescribed strain changes the affine offset, not dσ/du.
        // The reference transpose contains the same original tensor and B-bar
        // projection, without admitting the total thermal strain as elastic.
        let displacement = self.reference_bulk_stress_pullback(&scaled, &mut poll)?;
        checkpoint(&mut poll)?;
        Ok(BulkStressPullback3 {
            displacement,
            scales,
        })
    }
}
