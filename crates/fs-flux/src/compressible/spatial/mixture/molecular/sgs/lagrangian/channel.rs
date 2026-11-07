//! Physical-wall geometry for the shared Lagrangian mixture owner.

use super::*;
use crate::compressible::spatial::molecular::sgs::channel_axis;

impl SubgridTransport {
    pub(in crate::compressible::spatial::mixture::molecular) fn new_lagrangian_channel(
        flow: &MixtureEuler3d,
        config: MixtureLagrangianSubgridConfig,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries).map_err(dynamic::dynamic_error)?;
        Self::new_lagrangian_owner(flow, LagrangianConfig::FixedTransport(config), Some(normal))
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn new_lagrangian_channel_species(
        flow: &MixtureEuler3d,
        config: MixtureLagrangianSpeciesSubgridConfig,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries).map_err(dynamic::dynamic_error)?;
        Self::new_lagrangian_owner(flow, LagrangianConfig::Species(config), Some(normal))
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_channel_normal_axis(
        &self,
    ) -> Option<usize> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.normal_axis,
            _ => None,
        }
    }
}

impl LagrangianTransport {
    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn admit_wall(
        &self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<(), DissipationError> {
        if self.normal_axis == Some(side.axis()) && !matches!(wall, MolecularWall3d::NoSlip { .. })
        {
            return Err(DissipationError::Invalid(
                "Lagrangian channel requires two explicit no-slip walls",
            ));
        }
        Ok(())
    }
}
