//! Point response on the original exact-volume tangent space.

use super::{
    Cx, MixedImpulseControls, MixedViscoConstrainedImpulse, MixedViscoTetError,
    MixedViscoTetMotion, SolidPointMobility, SolidPointMobilityMatrix, Tet10ImpulsePoint,
    VolumetricTreatment, checkpoint, inertia_error, tet10,
};
use crate::hyper3::mixed::impulse::ConstrainedImpulseResponse;
use fs_sparse::Inertia;

/// One borrowed exact-volume endpoint and its reusable native KKT factor.
///
/// Pressure participates in every response column and the eventual impulse.
/// Mobility neither copies Maxwell histories nor subtracts a background
/// velocity. The borrow keeps geometry, mass, supports and accepted histories
/// fixed for the lifetime of the factor; a successor needs its own preparation.
#[derive(Debug)]
pub struct MixedViscoConstrainedResponse<'motion, 'material> {
    motion: &'motion MixedViscoTetMotion<'material>,
    response: ConstrainedImpulseResponse<'motion>,
    controls: MixedImpulseControls,
}

impl<'material> MixedViscoTetMotion<'material> {
    /// Prepare contact response with the original weak incompressibility rows.
    ///
    /// The accepted velocity must already satisfy the weak-volume-rate
    /// allowance. Current cofactors, mass and pressure constraints are assembled
    /// and factored once, independently of the number of contact directions.
    /// Original physical pressure and every material history remain borrowed.
    /// Finite-compliance motion uses its existing mass-only point methods.
    pub fn constrained_impulse_response<'motion>(
        &'motion self,
        held: &'motion [usize],
        pressure_gauge_node: Option<usize>,
        controls: MixedImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoConstrainedResponse<'motion, 'material>, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric != VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::FiniteCompliancePressureImpulse);
        }
        let response = ConstrainedImpulseResponse::prepare(
            &self.solid.mesh,
            &self.mass,
            held,
            pressure_gauge_node,
            controls,
            |e, q| {
                let f = &self.solid.states[e * self.solid.mesh.points_per_element + q]
                    .deformation_gradient;
                // Original endpoint derivative of det(F); no material trial,
                // Hessian evaluation or relaxation clock is introduced.
                Ok(core::array::from_fn(|ij| {
                    let (i, j) = (ij / 3, ij % 3);
                    let (a, b, c, d) = ((i + 1) % 3, (i + 2) % 3, (j + 1) % 3, (j + 2) % 3);
                    f[3 * a + c] * f[3 * b + d] - f[3 * a + d] * f[3 * b + c]
                }))
            },
            cx,
        )?;
        response.validate_velocity(&self.velocity_mps, cx)?;
        Ok(MixedViscoConstrainedResponse {
            motion: self,
            response,
            controls,
        })
    }
}

impl<'material> MixedViscoConstrainedResponse<'_, 'material> {
    /// Actual inertia of the one retained exact-constraint factor.
    pub fn kkt_inertia(&self) -> Inertia {
        self.response.kkt_inertia()
    }

    /// Complete pressure-constrained, cross-point/direction mobility [1/kg].
    ///
    /// Unit world directions use all ten signed P2 weights. Each column solves
    /// the homogeneous volume tangent equation directly through the same factor.
    /// One mesh-sized response is retained at a time. Original row/workspace
    /// bounds apply before allocation; cross terms are not repaired or clipped.
    pub fn point_mobility_matrix_per_kg(
        &self,
        rows: &[(Tet10ImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, MixedViscoTetError> {
        tet10::mobility_matrix_with(
            self.motion.velocity_mps.len(),
            self.motion.solid.elements,
            rows,
            maximum_rows,
            self.controls.balance,
            |loads| {
                self.response
                    .free_velocity_response(loads, cx)
                    .map_err(MixedViscoTetError::from)
            },
            |error| MixedViscoTetError::from(inertia_error(error)),
            cx,
        )
    }

    /// Directional point mobility with the same exact pressure constraint.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        point: Tet10ImpulsePoint,
        direction: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, MixedViscoTetError> {
        let mut matrix = self.point_mobility_matrix_per_kg(&[(point, direction)], 1, cx)?;
        Ok(SolidPointMobility {
            inverse_mass_per_kg: matrix.matrix_per_kg[0],
            solve: matrix.solves.remove(0),
        })
    }

    /// Sum original point loads before one exact-volume native event.
    /// Pressure impulse [Pa s], pressure zero-work and support work remain in
    /// the original report; retained physical pressure [Pa] is never replaced.
    pub fn point_impulses_candidate(
        &self,
        impulses: &[(Tet10ImpulsePoint, [f64; 3])],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<MixedViscoConstrainedImpulse<'material>, MixedViscoTetError> {
        let entries = tet10::distributions(
            self.motion.solid.elements,
            self.motion.velocity_mps.len(),
            impulses,
            self.controls.balance,
            cx,
        )
        .map_err(inertia_error)?;
        self.nodal_impulse_candidate(&entries, origin_m, cx)
    }

    /// One point impulse using the same signed trace and retained factor.
    pub fn point_impulse_candidate(
        &self,
        point: Tet10ImpulsePoint,
        impulse_world: [f64; 3],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<MixedViscoConstrainedImpulse<'material>, MixedViscoTetError> {
        let entries = tet10::distribution(
            point,
            self.motion.solid.elements,
            self.motion.velocity_mps.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.nodal_impulse_candidate(&entries, origin_m, cx)
    }

    pub(super) fn nodal_impulse_candidate(
        &self,
        loads: &[(usize, f64)],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<MixedViscoConstrainedImpulse<'material>, MixedViscoTetError> {
        let positions = self.motion.impulse_positions(cx)?;
        let (velocity_mps, impulse) = self.response.impulse_candidate(
            &positions,
            &self.motion.velocity_mps,
            loads,
            origin_m,
            cx,
        )?;
        Ok(MixedViscoConstrainedImpulse {
            motion: self.motion.with_impulse_velocity(velocity_mps, cx)?,
            impulse,
        })
    }
}
