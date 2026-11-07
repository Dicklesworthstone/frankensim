//! The same smooth-wall SGS traction on the original species/energy owner.

use super::*;

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Add explicit equilibrium smooth-wall stress to one already configured
    /// adiabatic no-slip wall. Original molecular source coefficients and
    /// physical wall velocity are resampled on every stage and retry.
    pub fn with_equilibrium_sgs_wall(
        mut self,
        side: BoundarySide,
        config: EquilibriumWallConfig,
    ) -> Result<Self, DissipationError> {
        let index = side_index(side);
        config
            .admit_wall(self.walls[index])
            .map_err(DissipationError::Invalid)?;
        let sgs = self.sgs.as_mut().ok_or(DissipationError::Invalid(
            "equilibrium wall stress requires an explicit SGS operator",
        ))?;
        sgs.wall_models[index] = Some(config);
        Ok(self)
    }

    /// Exact wall controls, independent of local/dynamic momentum or scalars.
    pub fn equilibrium_sgs_wall(&self, side: BoundarySide) -> Option<EquilibriumWallConfig> {
        self.sgs.as_ref()?.wall_models[side_index(side)]
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::mixture::molecular) fn has_equilibrium_wall(
        &self,
    ) -> bool {
        self.wall_models.iter().any(Option::is_some)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::mixture::molecular) fn add_wall(
        &self,
        side: BoundarySide,
        wall: MolecularWall3d,
        primitive: Primitive,
        distance: f64,
        species: usize,
        mut receipt: TransportFaceReport,
        out: &mut [f64],
    ) -> Result<TransportFaceReport, DissipationError> {
        let Some(config) = self.wall_models[side_index(side)] else {
            return Ok(receipt);
        };
        let molecular = receipt.dynamic_viscosity_pa_s;
        let viscosity = config
            .viscosity(
                wall,
                side.axis(),
                primitive.velocity,
                primitive.density,
                primitive.sound_speed,
                distance,
                molecular,
            )
            .map_err(DissipationError::Invalid)?;
        let MolecularWall3d::NoSlip { velocity, .. } = wall else {
            unreachable!();
        };
        // The source has already supplied impermeability, normal stress and
        // molecular wall work. Add only the missing tangential traction and
        // its work at the SAME physical wall velocity, in normal/tangent order.
        let extra = viscosity.secant - molecular;
        for (component, axis) in tangents(side.axis()).into_iter().enumerate() {
            let momentum =
                extra * side.sign() * (primitive.velocity[axis] - velocity[axis]) / distance;
            out[species + 1 + component] += momentum;
            out[species + 3] += momentum * velocity[axis];
        }
        // Preserve the source bound and add the nonlinear wall Jacobian's
        // excess over its molecular shear contribution. The owner subsequently
        // supplies the half-cell and three-axis inverse-length factors.
        receipt.diffusivity_m2_s += 4.0 * (viscosity.tangent - molecular) / primitive.density;
        if !out.iter().all(|v| v.is_finite()) || !receipt.diffusivity_m2_s.is_finite() {
            return Err(DissipationError::Invalid(
                "unrepresentable equilibrium SGS wall flux",
            ));
        }
        Ok(receipt)
    }
}
